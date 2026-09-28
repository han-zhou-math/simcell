! Stage 12 persistent one-way FSI -> actin manager.
!
! This module adds no Cartesian PDE discretization and no C/C++ entry point.
! A step is deliberately only
!
!   accepted FSI snapshot
!     -> Stage 10 snapshot-to-actin coefficients
!     -> Stage 11 packed boundary GMRES
!     -> Stage 10 scalar extension
!     -> unchanged Stage 06 semi-periodic C++ multigrid backend.
!
! The Stage-12 manager itself remains a one-way transport consumer.  Stage 14
! reuses it without changing the actin transaction and reads only its accepted
! state through get_actin_oneway_feedback_state before the next FSI solve.
! The manager's purpose is lifetime and transaction control.  It owns accepted
! theta_n/theta_c histories and consumes each monotonically numbered
! velocity/geometry snapshot at most once.
module fsi_actin_oneway_harness_mod
  use parameters, only: dp,nx,ny,npts,xmin,xmax,ymin,ymax,hg,cpi,pi,tupi,scc, &
       stage12_actin_dw,stage12_actin_theta0, &
       stage12_initial_network_concentration,stage12_initial_free_concentration, &
       enforce_fixed_timestep,nfreq, &
       stage12_localized_polymerization,stage12_pnas_balanced_actin_profile
  use grid_types, only: LagrangianGrid,EulerianGrid
  use geometry_mod, only: update_geometry_wrapper
  use fsi_transport_snapshot_mod, only: fsi_transport_snapshot_t, &
       FSI_SNAPSHOT_OK,copy_snapshot_marker_geometry
  use actin_model_mod, only: actin_coefficients_t
  use scalar_operator_mod, only: is_finite_scalar
  use generic_gmres_mod, only: gmres_options_t
  use fsi_actin_sequential_harness_mod, only: ACTIN_SEQUENTIAL_OK, &
       actin_species_state_t,actin_sequential_state_t, &
       actin_sequential_problem_t,actin_real_evaluator_context_t, &
       prepare_actin_problem_from_snapshot,evaluate_actin_species_real
  use fsi_actin_block_harness_mod, only: ACTIN_BLOCK_OK, &
       actin_block_result_t,advance_actin_block
  implicit none
  private

  integer,parameter,public::ACTIN_ONEWAY_OK=0
  integer,parameter,public::ACTIN_ONEWAY_INVALID=1
  integer,parameter,public::ACTIN_ONEWAY_NOT_READY=2
  integer,parameter,public::ACTIN_ONEWAY_STALE=3
  integer,parameter,public::ACTIN_ONEWAY_SOLVE_FAILED=4

  type,public::actin_oneway_manager_t
    private
    logical::initialized=.false.
    integer::accepted_snapshot_id=0
    type(actin_coefficients_t)::coefficients
    type(actin_sequential_state_t)::state
    type(LagrangianGrid)::lag_geometry
    type(EulerianGrid)::tagged_geometry
  end type actin_oneway_manager_t

  public::initialize_actin_oneway,advance_actin_oneway
  public::get_actin_oneway_mass,get_actin_oneway_fields
  public::get_actin_oneway_snapshot_id
  public::get_actin_oneway_feedback_state

contains

  subroutine initialize_actin_oneway(manager,marker_x,marker_y,coeff,status)
    type(actin_oneway_manager_t),intent(out)::manager
    real(dp),intent(in)::marker_x(:),marker_y(:)
    type(actin_coefficients_t),intent(in)::coeff
    integer,intent(out)::status
    integer::marker
    integer::nearest(nx,ny)

    status=ACTIN_ONEWAY_INVALID
    if(size(marker_x)/=npts .or. size(marker_y)/=npts)return
    if(.not.all(is_finite_scalar(marker_x)) .or. &
         .not.all(is_finite_scalar(marker_y)))return
    ! The correction and trace stencils require the material interface to be
    ! strictly separated from the outer box.  Stage 12 does not add contact.
    if(any(marker_x<=xmin+2.0_dp*hg) .or. &
         any(marker_x>=xmax-2.0_dp*hg))return
    if(any(marker_y<=ymin+2.0_dp*hg) .or. &
         any(marker_y>=ymax-2.0_dp*hg))return
    if(.not.valid_coefficients(coeff))return

    pi=cpi
    tupi=2.0_dp*cpi
    do marker=1,npts
      scc(marker)=tupi*real(marker-1,dp)/real(npts,dp)
    enddo
    call manager%lag_geometry%init(npts)
    call manager%tagged_geometry%init(nx,ny,xmin,xmax,ymin,ymax)
    manager%lag_geometry%x=marker_x
    manager%lag_geometry%y=marker_y
    manager%lag_geometry%x_old=marker_x
    manager%lag_geometry%y_old=marker_y
    call update_geometry_wrapper(manager%lag_geometry, &
         manager%tagged_geometry,1.0_dp,-1)
    call compute_full_nearest_map(marker_x,marker_y, &
         manager%tagged_geometry,nearest)
    ! Legacy tag_grid stores only a near-interface map.  The moving-history
    ! source needs the accepted full map, so the persistent owner supplies it.
    manager%tagged_geometry%kaio(1:nx,1:ny)=nearest

    call initialize_uniform_species(manager%state%network,stage12_initial_network_concentration)
    call initialize_uniform_species(manager%state%free,stage12_initial_free_concentration)
    allocate(manager%state%geometry%marker_x(npts), &
         manager%state%geometry%marker_y(npts), &
         manager%state%geometry%nearest_marker(nx,ny))
    manager%state%geometry%marker_x=marker_x
    manager%state%geometry%marker_y=marker_y
    manager%state%geometry%nearest_marker=nearest
    manager%coefficients=coeff
    manager%accepted_snapshot_id=0
    manager%initialized=.true.
    status=ACTIN_ONEWAY_OK
  end subroutine initialize_actin_oneway

  subroutine advance_actin_oneway(manager,snapshot,snapshot_id,dt,status)
    type(actin_oneway_manager_t),intent(inout)::manager
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    integer,intent(in)::snapshot_id
    real(dp),intent(in)::dt
    integer,intent(out)::status
    type(LagrangianGrid)::candidate_lag
    type(EulerianGrid)::candidate_eul
    type(actin_sequential_state_t)::candidate_state
    type(actin_sequential_problem_t)::problem
    type(actin_real_evaluator_context_t)::context
    type(actin_block_result_t)::block_result
    real(dp)::marker_x(npts),marker_y(npts),snapshot_normal(npts,2)
    integer::nearest(nx,ny),snapshot_status,solve_status

    status=ACTIN_ONEWAY_NOT_READY
    if(.not.manager%initialized)return
    status=ACTIN_ONEWAY_INVALID
    if(.not.is_finite_scalar(dt) .or. dt<=0.0_dp)return
    ! Exact monotone ownership prevents an old velocity/geometry value object
    ! from being silently applied twice to a newer accepted actin history.
    if(snapshot_id/=manager%accepted_snapshot_id+1)then
      status=ACTIN_ONEWAY_STALE
      return
    endif
    call copy_snapshot_marker_geometry(snapshot,marker_x,marker_y, &
         snapshot_normal,snapshot_id,snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)return
    if(any(marker_x<=xmin+2.0_dp*hg) .or. &
         any(marker_x>=xmax-2.0_dp*hg))return
    if(any(marker_y<=ymin+2.0_dp*hg) .or. &
         any(marker_y>=ymax-2.0_dp*hg))return

    ! Work only on disposable deep copies until the packed solve and its
    ! explicit final pair have passed.  Intrinsic assignment deep-copies all
    ! allocatable grid and species components.
    candidate_lag=manager%lag_geometry
    candidate_eul=manager%tagged_geometry
    candidate_state=manager%state
    candidate_lag%x_old=manager%state%geometry%marker_x
    candidate_lag%y_old=manager%state%geometry%marker_y
    candidate_lag%x=marker_x
    candidate_lag%y=marker_y
    call update_geometry_wrapper(candidate_lag,candidate_eul,dt,0)
    ! The snapshot is the single normal/velocity authority for this step.
    ! Reusing its normals also makes the Robin beta and trace reconstruction
    ! refer to the same orientation.
    candidate_lag%normal=snapshot_normal
    call compute_full_nearest_map(marker_x,marker_y,candidate_eul,nearest)
    ! tag_grid overwrites kaio with the current map.  Restore the accepted map
    ! so C^m is evaluated at X^{m+1} about its old polynomial center X^m on
    ! every swept cell:
    !       C^m(X^{m+1}; center=X^m)/dt.
    candidate_eul%kaio(1:nx,1:ny)= &
         manager%state%geometry%nearest_marker

    call prepare_actin_problem_from_snapshot(snapshot,snapshot_id, &
         manager%coefficients,dt,manager%state%free%traces%value_average, &
         nearest,problem,solve_status)
    if(solve_status/=ACTIN_SEQUENTIAL_OK)then
      status=ACTIN_ONEWAY_INVALID
      return
    endif
    problem%gmres_options=gmres_options_t(restart_length=40, &
         max_iterations=160,relative_tolerance=1.0e-9_dp, &
         absolute_tolerance=1.0e-10_dp)
    context%lag_geometry=candidate_lag
    context%tagged_geometry=candidate_eul
    context%evaluations=0
    call advance_actin_block(problem,candidate_state, &
         evaluate_actin_species_real,context,solve_status,block_result)
    if(solve_status/=ACTIN_BLOCK_OK)then
      status=ACTIN_ONEWAY_SOLVE_FAILED
      return
    endif

    ! The only accepted writes in the advance.  All failure paths above leave
    ! the prior two species, geometry, and snapshot id unchanged.
    manager%state=candidate_state
    manager%lag_geometry=candidate_lag
    manager%tagged_geometry=candidate_eul
    manager%accepted_snapshot_id=snapshot_id
    status=ACTIN_ONEWAY_OK
    call report_actin_gmres(snapshot_id,block_result)
    call report_figure2_actin_diagnostics(snapshot_id,marker_x,marker_y, &
         problem,manager%state,manager%coefficients)
  end subroutine advance_actin_oneway

  subroutine report_actin_gmres(snapshot_id,result)
    integer,intent(in)::snapshot_id
    type(actin_block_result_t),intent(in)::result
    integer::converged
    real(dp)::final_residual_max

    converged=merge(1,0,result%gmres%converged)
    final_residual_max=0.0_dp
    if(allocated(result%final_residual)) &
         final_residual_max=maxval(abs(result%final_residual))
    write(*,'(a,1x,6(i0,1x),4(es24.16,1x))') 'FIG2_ACTIN_GMRES', &
         snapshot_id,converged,result%status,result%evaluator_status, &
         result%gmres%iterations,result%gmres%restarts, &
         result%gmres%initial_residual_norm,result%gmres%absolute_residual, &
         result%gmres%relative_residual,final_residual_max
  end subroutine report_actin_gmres

  subroutine report_figure2_actin_diagnostics(snapshot_id,marker_x,marker_y, &
       problem,state,coeff)
    integer,intent(in)::snapshot_id
    real(dp),intent(in)::marker_x(:),marker_y(:)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(in)::state
    type(actin_coefficients_t),intent(in)::coeff
    real(dp)::dual_length(npts),profile(npts),old_g(npts),new_g(npts)
    real(dp)::rate_coefficient(npts),realized_flux(npts)
    real(dp)::normal_gradient(npts),relative_advective_normal(npts)
    real(dp)::network_relative_flux(npts),boundary_residual(npts)
    real(dp)::same_time_factor(npts),weight(npts)
    real(dp)::profile_integral,realized_integral,same_time_integral
    real(dp)::lagged_saturation,descriptive_saturation,weight_sum
    real(dp)::denominator,s
    integer::marker,previous_marker,next_marker,same_time_valid

    ! These records are emitted only after the complete packed actin solve has
    ! committed.  They use the same average traces as the affine Robin system,
    ! so `realized_flux` is the flux actually applied by the code rather than
    ! the continuum same-time approximation.
    do marker=1,npts
      previous_marker=modulo(marker-2,npts)+1
      next_marker=modulo(marker,npts)+1
      dual_length(marker)=0.5_dp*( &
           hypot(marker_x(marker)-marker_x(previous_marker), &
                 marker_y(marker)-marker_y(previous_marker))+ &
           hypot(marker_x(next_marker)-marker_x(marker), &
                 marker_y(next_marker)-marker_y(marker)))
      profile(marker)=1.0_dp
      if(stage12_localized_polymerization .and. &
           .not.stage12_pnas_balanced_actin_profile)then
        ! Match prepare_actin_problem_from_snapshot exactly: the localized
        ! profile is evaluated in material-marker coordinates.
        s=2.0_dp*cpi*real(marker-1,dp)/real(npts,dp)
        profile(marker)=2.0_dp-tanh(s**6/stage12_actin_dw)- &
             tanh((2.0_dp*cpi-s)**6/stage12_actin_dw)
      endif
    enddo

    old_g=problem%frozen_old_free_trace
    new_g=state%free%traces%value_average
    rate_coefficient=problem%coupling_q*coeff%network_diffusivity
    realized_flux=rate_coefficient*new_g
    normal_gradient=state%network%traces%normal_average+ &
         0.5_dp*state%network%density
    do marker=1,npts
      relative_advective_normal(marker)=dot_product( &
           problem%network_operator%velocity_marker(marker,:)- &
           problem%interface_velocity(marker,:), &
           problem%marker_normal(marker,:))
    enddo
    network_relative_flux=relative_advective_normal* &
         state%network%traces%value_average- &
         coeff%network_diffusivity*normal_gradient
    boundary_residual=network_relative_flux+realized_flux

    same_time_valid=1
    same_time_factor=0.0_dp
    do marker=1,npts
      denominator=new_g(marker)+stage12_actin_theta0
      if(abs(denominator)<=100.0_dp*epsilon(1.0_dp))then
        same_time_valid=0
      else
        same_time_factor(marker)=new_g(marker)/denominator
      endif
    enddo
    weight=profile*dual_length
    weight_sum=sum(weight)
    profile_integral=sum(weight)
    realized_integral=sum(realized_flux*dual_length)
    same_time_integral=0.0_dp
    lagged_saturation=0.0_dp
    descriptive_saturation=0.0_dp
    if(stage12_localized_polymerization .and. &
         .not.stage12_pnas_balanced_actin_profile .and. &
         weight_sum>100.0_dp*epsilon(1.0_dp))then
      if(same_time_valid==1)then
        same_time_integral=coeff%membrane_rate*sum(weight*same_time_factor)
        descriptive_saturation=sum(weight*same_time_factor)/weight_sum
      endif
      lagged_saturation=sum(weight*new_g/max(old_g+stage12_actin_theta0, &
           100.0_dp*epsilon(1.0_dp)))/weight_sum
    endif

    write(*,'(a,1x,2(i0,1x),7(es24.16,1x))') &
         'FIG2_ACTIN_POLYMERIZATION',snapshot_id,same_time_valid, &
         profile_integral,realized_integral,same_time_integral, &
         lagged_saturation,descriptive_saturation, &
         maxval(abs(boundary_residual)), &
         sqrt(sum(boundary_residual*boundary_residual))
    ! In fixed-step campaigns snapshot_id equals the committed outer step.
    if(enforce_fixed_timestep .and. nfreq>0)then
      if(mod(snapshot_id,nfreq)/=0)return
    endif
    do marker=1,npts
      write(*,'(a,1x,2(i0,1x),20(es24.16,1x))') 'FIG2_ACTIN_MARKER', &
           snapshot_id,marker,marker_x(marker),marker_y(marker), &
           dual_length(marker),profile(marker),old_g(marker),new_g(marker), &
           state%free%traces%value_interior(marker), &
           state%network%traces%value_average(marker), &
           state%network%traces%value_interior(marker), &
           normal_gradient(marker),rate_coefficient(marker), &
           realized_flux(marker),problem%marker_normal(marker,1), &
           problem%marker_normal(marker,2), &
           problem%network_operator%velocity_marker(marker,1), &
           problem%network_operator%velocity_marker(marker,2), &
           problem%interface_velocity(marker,1), &
           problem%interface_velocity(marker,2), &
           network_relative_flux(marker),boundary_residual(marker)
    enddo
  end subroutine report_figure2_actin_diagnostics

  subroutine get_actin_oneway_mass(manager,total_mass,status)
    type(actin_oneway_manager_t),intent(in)::manager
    real(dp),intent(out)::total_mass
    integer,intent(out)::status
    integer::i,j
    total_mass=0.0_dp
    status=ACTIN_ONEWAY_NOT_READY
    if(.not.manager%initialized)return
    do j=1,ny
      do i=1,nx
        if(manager%tagged_geometry%idf(i,j)/=1)cycle
        total_mass=total_mass+(manager%state%network%field(i,j)+ &
             manager%state%free%field(i,j))*hg*hg
      enddo
    enddo
    status=ACTIN_ONEWAY_OK
  end subroutine get_actin_oneway_mass

  subroutine get_actin_oneway_fields(manager,network,free,status)
    type(actin_oneway_manager_t),intent(in)::manager
    real(dp),intent(out)::network(:,:),free(:,:)
    integer,intent(out)::status
    integer::i,j
    network=0.0_dp; free=0.0_dp
    status=ACTIN_ONEWAY_NOT_READY
    if(.not.manager%initialized)return
    if(any(shape(network)/=[nx,ny]) .or. any(shape(free)/=[nx,ny]))then
      status=ACTIN_ONEWAY_INVALID
      return
    endif
    do j=1,ny
      do i=1,nx
        if(manager%tagged_geometry%idf(i,j)/=1)cycle
        network(i,j)=manager%state%network%field(i,j)
        free(i,j)=manager%state%free%field(i,j)
      enddo
    enddo
    status=ACTIN_ONEWAY_OK
  end subroutine get_actin_oneway_fields

  subroutine get_actin_oneway_snapshot_id(manager,snapshot_id,status)
    type(actin_oneway_manager_t),intent(in)::manager
    integer,intent(out)::snapshot_id,status
    snapshot_id=0
    status=ACTIN_ONEWAY_NOT_READY
    if(.not.manager%initialized)return
    snapshot_id=manager%accepted_snapshot_id
    status=ACTIN_ONEWAY_OK
  end subroutine get_actin_oneway_snapshot_id

  subroutine get_actin_oneway_feedback_state(manager,auxiliary_field, &
       correction_grid,cell_side,marker_trace,accepted_snapshot_id,status)
    type(actin_oneway_manager_t),intent(in)::manager
    real(dp),intent(out)::auxiliary_field(:,:),correction_grid(:,:)
    integer,intent(out)::cell_side(:,:)
    real(dp),intent(out)::marker_trace(:)
    integer,intent(out)::accepted_snapshot_id,status

    ! The getter returns value copies only.  In particular, no grid retagging,
    ! history update, or pointer alias can occur while the FSI outer iteration
    ! consumes this accepted time level.
    auxiliary_field=0.0_dp
    correction_grid=0.0_dp
    cell_side=0
    marker_trace=0.0_dp
    accepted_snapshot_id=0
    status=ACTIN_ONEWAY_NOT_READY
    if(.not.manager%initialized)return
    status=ACTIN_ONEWAY_INVALID
    if(any(shape(auxiliary_field)/=[nx,ny]))return
    if(any(shape(correction_grid)/=[nx,ny]))return
    if(any(shape(cell_side)/=[nx,ny]))return
    if(size(marker_trace)/=npts)return
    if(.not.allocated(manager%state%network%field))return
    if(.not.allocated(manager%state%network%correction_grid))return
    if(.not.allocated(manager%state%network%traces%value_interior))return
    if(any(shape(manager%state%network%field)/=[nx,ny]))return
    if(any(shape(manager%state%network%correction_grid)/=[nx,ny]))return
    if(size(manager%state%network%traces%value_interior)/=npts)return
    if(.not.all(is_finite_scalar(manager%state%network%field)))return
    if(.not.all(is_finite_scalar( &
         manager%state%network%correction_grid)))return
    if(.not.all(is_finite_scalar( &
         manager%state%network%traces%value_interior)))return
    if(any(abs(manager%tagged_geometry%idf(1:nx,1:ny))/=1))return

    auxiliary_field=manager%state%network%field
    correction_grid=manager%state%network%correction_grid
    cell_side=manager%tagged_geometry%idf(1:nx,1:ny)
    marker_trace=manager%state%network%traces%value_interior
    accepted_snapshot_id=manager%accepted_snapshot_id
    status=ACTIN_ONEWAY_OK
  end subroutine get_actin_oneway_feedback_state

  subroutine initialize_uniform_species(species,value)
    type(actin_species_state_t),intent(out)::species
    real(dp),intent(in)::value
    allocate(species%field(nx,ny),species%density(npts), &
         species%correction_coefficients(npts,6), &
         species%correction_grid(nx,ny), &
         species%traces%value_average(npts), &
         species%traces%normal_average(npts), &
         species%traces%value_interior(npts), &
         species%traces%value_exterior(npts), &
         species%traces%normal_interior(npts), &
         species%traces%normal_exterior(npts))
    species%field=value
    species%density=0.0_dp
    species%correction_coefficients=0.0_dp
    species%correction_grid=0.0_dp
    species%traces%value_average=value
    species%traces%normal_average=0.0_dp
    species%traces%value_interior=value
    species%traces%value_exterior=value
    species%traces%normal_interior=0.0_dp
    species%traces%normal_exterior=0.0_dp
  end subroutine initialize_uniform_species

  subroutine compute_full_nearest_map(marker_x,marker_y,eul,nearest)
    real(dp),intent(in)::marker_x(:),marker_y(:)
    type(EulerianGrid),intent(in)::eul
    integer,intent(out)::nearest(:,:)
    real(dp)::xc,yc,distance,best
    integer::i,j,marker
    do j=1,ny
      yc=eul%y_min+(real(j,dp)-0.5_dp)*eul%dy
      do i=1,nx
        xc=eul%x_min+(real(i,dp)-0.5_dp)*eul%dx
        best=huge(1.0_dp); nearest(i,j)=1
        do marker=1,npts
          distance=(xc-marker_x(marker))**2+(yc-marker_y(marker))**2
          if(distance<best)then
            best=distance
            nearest(i,j)=marker
          endif
        enddo
      enddo
    enddo
  end subroutine compute_full_nearest_map

  logical function valid_coefficients(coeff)
    type(actin_coefficients_t),intent(in)::coeff
    valid_coefficients=.false.
    if(.not.is_finite_scalar(coeff%network_transport_fraction))return
    if(.not.is_finite_scalar(coeff%network_diffusivity))return
    if(.not.is_finite_scalar(coeff%free_diffusivity))return
    if(.not.is_finite_scalar(coeff%turnover_rate))return
    if(.not.is_finite_scalar(coeff%membrane_rate))return
    if(coeff%network_transport_fraction<=0.0_dp .or. &
         coeff%network_transport_fraction>=1.0_dp)return
    if(coeff%network_diffusivity<=0.0_dp .or. &
         coeff%free_diffusivity<=0.0_dp)return
    if(coeff%turnover_rate<0.0_dp .or. coeff%membrane_rate<0.0_dp)return
    valid_coefficients=.true.
  end function valid_coefficients

end module fsi_actin_oneway_harness_mod
