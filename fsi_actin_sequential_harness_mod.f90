! Stage 10 network-then-free actin orchestration.
!
! Provenance: the callback/state separation was studied in code/dualchem at
! commit 29224de, but this routine is reimplemented for the required opposite
! order (network then free) and deliberately omits that source's block solver.
!
! This file is intentionally a harness, not a Cartesian PDE backend.  Both
! species are evaluated through actin_extension_mod, which normalizes by the
! physical diffusion and calls the unchanged Stage 06 semi-periodic C++ MG
! symbol.  This module owns only the boundary-density sequence
!
!   network/F-actin at m+1  ->  free/G-actin at m+1,
!
! and commits both trials atomically.  It contains no block residual and no
! actin-to-FSI feedback; those belong to later stages.
module fsi_actin_sequential_harness_mod
  use parameters, only: dp,nx,ny,npts,cpi,stage12_actin_dw, &
       stage12_actin_theta0,stage12_localized_polymerization, &
       stage12_pnas_balanced_actin_profile
  use fsi_transport_snapshot_mod, only: fsi_transport_snapshot_t, &
       FSI_SNAPSHOT_OK,copy_snapshot_face_velocity, &
       copy_snapshot_marker_velocity,copy_snapshot_interface_velocity, &
       copy_snapshot_marker_geometry
  use actin_model_mod, only: actin_coefficients_t,network_robin_beta
  use actin_pnas_profile_mod, only: pnas_actin_gaussian_profile
  use scalar_operator_mod, only: scalar_operator_t,scalar_rhs_t, &
       SCALAR_STATUS_OK,validate_scalar_operator,is_finite_scalar
  use generic_gmres_mod, only: gmres_options_t,gmres_result_t, &
       GMRES_STATUS_SUCCESS,solve_affine_gmres
  use actin_extension_mod, only: interface_traces_t,evaluate_actin_extension, &
       ACTIN_EXTENSION_OK
  use grid_types, only: LagrangianGrid,EulerianGrid
  implicit none
  private

  integer,parameter,public::ACTIN_SPECIES_NETWORK=1
  integer,parameter,public::ACTIN_SPECIES_FREE=2

  integer,parameter,public::ACTIN_SEQUENTIAL_OK=0
  integer,parameter,public::ACTIN_SEQUENTIAL_INVALID=1
  integer,parameter,public::ACTIN_SEQUENTIAL_NETWORK_SOLVE_FAILED=2
  integer,parameter,public::ACTIN_SEQUENTIAL_NETWORK_FINAL_FAILED=3
  integer,parameter,public::ACTIN_SEQUENTIAL_FREE_SOLVE_FAILED=4
  integer,parameter,public::ACTIN_SEQUENTIAL_FREE_FINAL_FAILED=5

  type,public::actin_extension_trial_t
    real(dp),allocatable::field(:,:)
    real(dp),allocatable::correction_coefficients(:,:)
    real(dp),allocatable::correction_grid(:,:)
    type(interface_traces_t)::traces
  end type actin_extension_trial_t

  type,public::actin_species_state_t
    real(dp),allocatable::field(:,:),density(:)
    real(dp),allocatable::correction_coefficients(:,:)
    real(dp),allocatable::correction_grid(:,:)
    type(interface_traces_t)::traces
  end type actin_species_state_t

  type,public::actin_geometry_history_t
    real(dp),allocatable::marker_x(:),marker_y(:)
    integer,allocatable::nearest_marker(:,:)
  end type actin_geometry_history_t

  type,public::actin_sequential_state_t
    type(actin_species_state_t)::network
    type(actin_species_state_t)::free
    type(actin_geometry_history_t)::geometry
  end type actin_sequential_state_t

  type,public::actin_sequential_problem_t
    type(actin_coefficients_t)::coefficients
    type(scalar_operator_t)::network_operator,free_operator
    real(dp),allocatable::beta_network(:),beta_free(:)
    real(dp),allocatable::frozen_old_free_trace(:)
    real(dp),allocatable::marker_normal(:,:)
    real(dp),allocatable::interface_velocity(:,:)
    real(dp),allocatable::coupling_q(:)
    real(dp)::dt=0.0_dp
    integer::snapshot_id=0
    type(gmres_options_t)::gmres_options
    type(actin_geometry_history_t)::next_geometry
  end type actin_sequential_problem_t

  type,public::actin_sequential_result_t
    integer::status=ACTIN_SEQUENTIAL_INVALID
    integer::evaluator_status=0
    type(gmres_result_t)::network_gmres,free_gmres
    real(dp),allocatable::network_final_residual(:),free_final_residual(:)
  end type actin_sequential_result_t

  ! The production evaluator binds disposable scalar solves to one already
  ! tagged geometry.  The orchestration itself remains testable with an affine
  ! evaluator and therefore does not hide ordering behind the C++ backend.
  type,public::actin_real_evaluator_context_t
    type(LagrangianGrid)::lag_geometry
    type(EulerianGrid)::tagged_geometry
    integer::evaluations=0
  end type actin_real_evaluator_context_t

  abstract interface
    subroutine actin_species_evaluator(species,density,operator,rhs,dt, &
         u_old,Cold,trial,evaluator_context,status)
      import::dp,scalar_operator_t,scalar_rhs_t,actin_extension_trial_t
      integer,intent(in)::species
      real(dp),intent(in)::density(:),dt,u_old(:,:),Cold(:,:)
      type(scalar_operator_t),intent(in)::operator
      type(scalar_rhs_t),intent(in)::rhs
      type(actin_extension_trial_t),intent(out)::trial
      class(*),intent(inout)::evaluator_context
      integer,intent(out)::status
    end subroutine actin_species_evaluator
  end interface

  public::actin_species_evaluator
  public::prepare_actin_problem_from_snapshot
  public::advance_actin_sequential
  public::evaluate_actin_species_real

contains

  subroutine prepare_actin_problem_from_snapshot(snapshot,snapshot_id,coeff, &
       dt,frozen_old_free_trace,next_nearest_marker,problem,status)
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    integer,intent(in)::snapshot_id
    type(actin_coefficients_t),intent(in)::coeff
    real(dp),intent(in)::dt,frozen_old_free_trace(:)
    integer,intent(in)::next_nearest_marker(:,:)
    type(actin_sequential_problem_t),intent(out)::problem
    integer,intent(out)::status
    real(dp)::x_new(npts),y_new(npts),normal(npts,2)
    real(dp)::free_marker_velocity(npts,2),interface_velocity(npts,2)
    real(dp)::vc_dot_n,interface_dot_n,a,s,profile,local_rate
    real(dp)::profile_values(npts),marker_metric(npts),mean_metric
    integer::marker,previous_marker,next_marker,snapshot_status,allocation_status

    status=ACTIN_SEQUENTIAL_INVALID
    if(snapshot_id<=0)return
    if(.not.is_finite_scalar(dt))return
    if(dt<=0.0_dp)return
    if(size(frozen_old_free_trace)/=npts)return
    if(any(shape(next_nearest_marker)/=[nx,ny]))return
    if(any(next_nearest_marker<1) .or. any(next_nearest_marker>npts))return
    if(.not.all(is_finite_scalar(frozen_old_free_trace)))return
    if(.not.valid_coefficients(coeff))return

    allocate(problem%network_operator%vx_face(nx+1,ny), &
         problem%network_operator%vy_face(nx,ny+1), &
         problem%network_operator%velocity_marker(npts,2), &
         problem%network_operator%flux_w(ny), &
         problem%network_operator%flux_e(ny), &
         problem%network_operator%flux_s(nx), &
         problem%network_operator%flux_n(nx), &
         problem%free_operator%vx_face(nx+1,ny), &
         problem%free_operator%vy_face(nx,ny+1), &
         problem%free_operator%velocity_marker(npts,2), &
         problem%free_operator%flux_w(ny),problem%free_operator%flux_e(ny), &
         problem%free_operator%flux_s(nx),problem%free_operator%flux_n(nx), &
         problem%beta_network(npts),problem%beta_free(npts), &
         problem%frozen_old_free_trace(npts),problem%marker_normal(npts,2), &
         problem%interface_velocity(npts,2),problem%coupling_q(npts), &
         problem%next_geometry%marker_x(npts), &
         problem%next_geometry%marker_y(npts), &
         problem%next_geometry%nearest_marker(nx,ny),stat=allocation_status)
    if(allocation_status/=0)return

    ! All face, marker, interface, and normal data below are copied from this
    ! one snapshot id.  The network velocity scales only v_c by a; V_Gamma is
    ! copied unchanged and appears solely in the Robin coefficients.
    call copy_snapshot_face_velocity(snapshot,problem%free_operator%vx_face, &
         problem%free_operator%vy_face,snapshot_id,snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)return
    call copy_snapshot_marker_geometry(snapshot,x_new,y_new,normal, &
         snapshot_id,snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)return
    call copy_snapshot_marker_velocity(snapshot,x_new,y_new, &
         free_marker_velocity,1.0_dp,snapshot_id,snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)return
    call copy_snapshot_interface_velocity(snapshot,interface_velocity, &
         snapshot_id,snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)return

    a=coeff%network_transport_fraction
    problem%free_operator%velocity_marker=free_marker_velocity
    problem%network_operator%vx_face=a*problem%free_operator%vx_face
    problem%network_operator%vy_face=a*problem%free_operator%vy_face
    problem%network_operator%velocity_marker=a*free_marker_velocity
    problem%network_operator%diffusion=coeff%network_diffusivity
    problem%network_operator%reaction=coeff%turnover_rate
    problem%free_operator%diffusion=coeff%free_diffusivity
    problem%free_operator%reaction=0.0_dp
    problem%network_operator%flux_w=0.0_dp
    problem%network_operator%flux_e=0.0_dp
    problem%network_operator%flux_s=0.0_dp
    problem%network_operator%flux_n=0.0_dp
    problem%free_operator%flux_w=0.0_dp
    problem%free_operator%flux_e=0.0_dp
    problem%free_operator%flux_s=0.0_dp
    problem%free_operator%flux_n=0.0_dp

    profile_values=0.0_dp
    marker_metric=1.0_dp
    if(stage12_pnas_balanced_actin_profile)then
      do marker=1,npts
        previous_marker=modulo(marker-2,npts)+1
        next_marker=modulo(marker,npts)+1
        marker_metric(marker)=0.5_dp*( &
             hypot(x_new(marker)-x_new(previous_marker), &
                   y_new(marker)-y_new(previous_marker))+ &
             hypot(x_new(next_marker)-x_new(marker), &
                   y_new(next_marker)-y_new(marker)))
        s=2.0_dp*cpi*real(marker-1,dp)/real(npts,dp)
        profile_values(marker)=pnas_actin_gaussian_profile(s,1.0_dp)
      enddo
      if(.not.all(is_finite_scalar(marker_metric)))return
      if(any(marker_metric<=100.0_dp*epsilon(1.0_dp)))return
      ! Center the material-coordinate profile before applying the PNAS
      ! |dX0/ds|/|dX/ds| factor.  This makes the discrete physical-arclength
      ! integral of the signed flux zero on every accepted interface.
      profile_values=profile_values-sum(profile_values)/real(npts,dp)
      mean_metric=sum(marker_metric)/real(npts,dp)
      profile_values=profile_values*mean_metric/marker_metric
      if(.not.all(is_finite_scalar(profile_values)))return
    endif

    do marker=1,npts
      vc_dot_n=dot_product(free_marker_velocity(marker,:),normal(marker,:))
      interface_dot_n=dot_product(interface_velocity(marker,:),normal(marker,:))
      local_rate=coeff%membrane_rate
      if(stage12_pnas_balanced_actin_profile)then
        local_rate=coeff%membrane_rate*profile_values(marker)
      elseif(stage12_localized_polymerization)then
        s=2.0_dp*cpi*real(marker-1,dp)/real(npts,dp)
        profile=2.0_dp-tanh(s**6/stage12_actin_dw)- &
             tanh((2.0_dp*cpi-s)**6/stage12_actin_dw)
        ! Yao--Li use j(s)*theta_c/(theta_c+theta_0).  Freezing only the
        ! denominator at the accepted trace preserves an affine block solve.
        local_rate=coeff%membrane_rate*profile/ &
             max(frozen_old_free_trace(marker)+stage12_actin_theta0, &
                 100.0_dp*epsilon(1.0_dp))
      endif
      problem%beta_network(marker)=network_robin_beta(coeff,vc_dot_n, &
           interface_dot_n)
      problem%beta_free(marker)=(interface_dot_n-vc_dot_n+local_rate)/ &
           coeff%free_diffusivity
      problem%coupling_q(marker)=local_rate/coeff%network_diffusivity
    enddo
    if(.not.all(is_finite_scalar(problem%beta_network)))return
    if(.not.all(is_finite_scalar(problem%beta_free)))return

    problem%coefficients=coeff
    problem%dt=dt
    problem%snapshot_id=snapshot_id
    problem%frozen_old_free_trace=frozen_old_free_trace
    problem%marker_normal=normal
    problem%interface_velocity=interface_velocity
    problem%next_geometry%marker_x=x_new
    problem%next_geometry%marker_y=y_new
    problem%next_geometry%nearest_marker=next_nearest_marker
    status=ACTIN_SEQUENTIAL_OK
  end subroutine prepare_actin_problem_from_snapshot

  subroutine advance_actin_sequential(problem,state,evaluate_species, &
       evaluator_context,status,result)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(inout)::state
    procedure(actin_species_evaluator)::evaluate_species
    class(*),intent(inout)::evaluator_context
    integer,intent(out)::status
    type(actin_sequential_result_t),intent(out),optional::result

    type(actin_sequential_state_t)::staged_state
    type(actin_extension_trial_t)::network_trial,free_trial
    type(actin_sequential_result_t)::work_result
    real(dp),allocatable::network_density(:),free_density(:)
    integer::nmarker,evaluator_status,allocation_status

    call initialize_result(work_result)
    status=ACTIN_SEQUENTIAL_INVALID
    if(.not.valid_problem_state(problem,state))then
      if(present(result))result=work_result
      return
    endif
    nmarker=size(state%network%density)
    allocate(network_density(nmarker),free_density(nmarker), &
         work_result%network_final_residual(nmarker), &
         work_result%free_final_residual(nmarker),stat=allocation_status)
    if(allocation_status/=0)then
      if(present(result))result=work_result
      return
    endif
    network_density=0.0_dp; free_density=0.0_dp
    work_result%network_final_residual=0.0_dp
    work_result%free_final_residual=0.0_dp
    staged_state=state

    ! First solve: theta_n uses only accepted network history and the accepted
    ! old theta_c trace frozen on the new interface.  No free trial exists yet.
    call solve_affine_gmres(network_callback,evaluator_context, &
         network_density,problem%gmres_options,work_result%network_gmres)
    if(.not.work_result%network_gmres%converged .or. &
         work_result%network_gmres%status/=GMRES_STATUS_SUCCESS)then
      work_result%status=ACTIN_SEQUENTIAL_NETWORK_SOLVE_FAILED
      work_result%evaluator_status=work_result%network_gmres%callback_status
      status=work_result%status
      if(present(result))result=work_result
      return
    endif
    call evaluate_network_residual(problem,state,network_density, &
         evaluate_species,evaluator_context,work_result%network_final_residual, &
         network_trial,evaluator_status)
    if(evaluator_status/=0 .or. .not.final_residual_is_accepted( &
         work_result%network_final_residual,work_result%network_gmres, &
         problem%gmres_options))then
      work_result%status=ACTIN_SEQUENTIAL_NETWORK_FINAL_FAILED
      work_result%evaluator_status=evaluator_status
      status=work_result%status
      if(present(result))result=work_result
      return
    endif
    call stage_species(staged_state%network,network_density,network_trial)

    ! Second solve: theta_c uses gamma times the just-staged theta_n^{m+1}
    ! field.  The public state still holds time m until both solves finalize.
    call solve_affine_gmres(free_callback,evaluator_context,free_density, &
         problem%gmres_options,work_result%free_gmres)
    if(.not.work_result%free_gmres%converged .or. &
         work_result%free_gmres%status/=GMRES_STATUS_SUCCESS)then
      work_result%status=ACTIN_SEQUENTIAL_FREE_SOLVE_FAILED
      work_result%evaluator_status=work_result%free_gmres%callback_status
      status=work_result%status
      if(present(result))result=work_result
      return
    endif
    call evaluate_free_residual(problem,state,staged_state%network, &
         free_density,evaluate_species,evaluator_context, &
         work_result%free_final_residual,free_trial,evaluator_status)
    if(evaluator_status/=0 .or. .not.final_residual_is_accepted( &
         work_result%free_final_residual,work_result%free_gmres, &
         problem%gmres_options))then
      work_result%status=ACTIN_SEQUENTIAL_FREE_FINAL_FAILED
      work_result%evaluator_status=evaluator_status
      status=work_result%status
      if(present(result))result=work_result
      return
    endif
    call stage_species(staged_state%free,free_density,free_trial)
    staged_state%geometry=problem%next_geometry

    ! The only accepted-state write in the routine.  Every failure above
    ! returns with both species, all traces, and both correction histories at m.
    state=staged_state
    work_result%status=ACTIN_SEQUENTIAL_OK
    work_result%evaluator_status=0
    status=ACTIN_SEQUENTIAL_OK
    if(present(result))result=work_result

  contains
    subroutine network_callback(density,residual,context,callback_status)
      real(dp),intent(in)::density(:)
      real(dp),intent(out)::residual(:)
      class(*),intent(inout)::context
      integer,intent(out)::callback_status
      type(actin_extension_trial_t)::disposable
      call evaluate_network_residual(problem,state,density,evaluate_species, &
           context,residual,disposable,callback_status)
    end subroutine network_callback

    subroutine free_callback(density,residual,context,callback_status)
      real(dp),intent(in)::density(:)
      real(dp),intent(out)::residual(:)
      class(*),intent(inout)::context
      integer,intent(out)::callback_status
      type(actin_extension_trial_t)::disposable
      call evaluate_free_residual(problem,state,staged_state%network,density, &
           evaluate_species,context,residual,disposable,callback_status)
    end subroutine free_callback
  end subroutine advance_actin_sequential

  subroutine evaluate_network_residual(problem,state,density,evaluate_species, &
       context,residual,trial,status)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(in)::state
    real(dp),intent(in)::density(:)
    procedure(actin_species_evaluator)::evaluate_species
    class(*),intent(inout)::context
    real(dp),intent(out)::residual(:)
    type(actin_extension_trial_t),intent(out)::trial
    integer,intent(out)::status
    type(scalar_rhs_t)::rhs
    integer::i,nmarker,nx_local,ny_local,allocation_status
    real(dp)::term,value,next_term
    logical::ok

    residual=0.0_dp; status=ACTIN_SEQUENTIAL_INVALID
    nmarker=size(density); nx_local=size(state%network%field,1)
    ny_local=size(state%network%field,2)
    allocate(rhs%volume(nx_local,ny_local),rhs%source_jump(nmarker), &
         stat=allocation_status)
    if(allocation_status/=0)return
    rhs%volume=0.0_dp; rhs%source_jump=0.0_dp
    call evaluate_species(ACTIN_SPECIES_NETWORK,density, &
         problem%network_operator,rhs,problem%dt,state%network%field, &
         state%network%correction_coefficients,trial,context,status)
    if(status/=0)return
    if(.not.valid_trial(trial,nx_local,ny_local,nmarker))then
      status=ACTIN_SEQUENTIAL_INVALID; return
    endif
    do i=1,nmarker
      call checked_multiply(0.5_dp,density(i),term,ok); if(.not.ok)return
      call checked_add(term,trial%traces%normal_average(i),next_term,ok)
      if(.not.ok)return
      term=next_term
      call checked_multiply(problem%beta_network(i), &
           trial%traces%value_average(i),value,ok); if(.not.ok)return
      call checked_add(term,value,next_term,ok); if(.not.ok)return
      term=next_term
      call checked_multiply(problem%coupling_q(i), &
           problem%frozen_old_free_trace(i),value,ok); if(.not.ok)return
      call checked_add(term,-value,residual(i),ok); if(.not.ok)return
    enddo
    if(.not.all(is_finite_scalar(residual)))then
      residual=0.0_dp; status=ACTIN_SEQUENTIAL_INVALID; return
    endif
    status=0
  end subroutine evaluate_network_residual

  subroutine evaluate_free_residual(problem,accepted_state,staged_network, &
       density,evaluate_species,context,residual,trial,status)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(in)::accepted_state
    type(actin_species_state_t),intent(in)::staged_network
    real(dp),intent(in)::density(:)
    procedure(actin_species_evaluator)::evaluate_species
    class(*),intent(inout)::context
    real(dp),intent(out)::residual(:)
    type(actin_extension_trial_t),intent(out)::trial
    integer,intent(out)::status
    type(scalar_rhs_t)::rhs
    integer::i,j,nmarker,nx_local,ny_local,allocation_status
    real(dp)::term,value,next_term
    logical::ok

    residual=0.0_dp; status=ACTIN_SEQUENTIAL_INVALID
    nmarker=size(density); nx_local=size(accepted_state%free%field,1)
    ny_local=size(accepted_state%free%field,2)
    allocate(rhs%volume(nx_local,ny_local),rhs%source_jump(nmarker), &
         stat=allocation_status)
    if(allocation_status/=0)return
    rhs%source_jump=0.0_dp
    do j=1,ny_local
      do i=1,nx_local
        call checked_multiply(problem%coefficients%turnover_rate, &
             staged_network%field(i,j),rhs%volume(i,j),ok)
        if(.not.ok)return
      enddo
    enddo
    call evaluate_species(ACTIN_SPECIES_FREE,density,problem%free_operator, &
         rhs,problem%dt,accepted_state%free%field, &
         accepted_state%free%correction_coefficients,trial,context,status)
    if(status/=0)return
    if(.not.valid_trial(trial,nx_local,ny_local,nmarker))then
      status=ACTIN_SEQUENTIAL_INVALID; return
    endif
    do i=1,nmarker
      call checked_multiply(0.5_dp,density(i),term,ok); if(.not.ok)return
      call checked_add(term,trial%traces%normal_average(i),next_term,ok)
      if(.not.ok)return
      term=next_term
      call checked_multiply(problem%beta_free(i), &
           trial%traces%value_average(i),value,ok); if(.not.ok)return
      call checked_add(term,value,residual(i),ok); if(.not.ok)return
    enddo
    if(.not.all(is_finite_scalar(residual)))then
      residual=0.0_dp; status=ACTIN_SEQUENTIAL_INVALID; return
    endif
    status=0
  end subroutine evaluate_free_residual

  subroutine evaluate_actin_species_real(species,density,operator,rhs,dt, &
       u_old,Cold,trial,evaluator_context,status)
    integer,intent(in)::species
    real(dp),intent(in)::density(:),dt,u_old(:,:),Cold(:,:)
    type(scalar_operator_t),intent(in)::operator
    type(scalar_rhs_t),intent(in)::rhs
    type(actin_extension_trial_t),intent(out)::trial
    class(*),intent(inout)::evaluator_context
    integer,intent(out)::status
    integer::backend_iterations

    status=ACTIN_SEQUENTIAL_INVALID
    if(species/=ACTIN_SPECIES_NETWORK .and. species/=ACTIN_SPECIES_FREE)return
    select type(context=>evaluator_context)
    type is(actin_real_evaluator_context_t)
      call evaluate_actin_extension(context%lag_geometry, &
           context%tagged_geometry,density,operator,rhs,dt,u_old,Cold, &
           trial%field,trial%correction_coefficients,trial%correction_grid, &
           trial%traces,status,backend_iterations)
      context%evaluations=context%evaluations+1
      if(status==ACTIN_EXTENSION_OK)status=0
    class default
      return
    end select
  end subroutine evaluate_actin_species_real

  subroutine stage_species(species,density,trial)
    type(actin_species_state_t),intent(inout)::species
    real(dp),intent(in)::density(:)
    type(actin_extension_trial_t),intent(in)::trial
    species%field=trial%field
    species%density=density
    species%correction_coefficients=trial%correction_coefficients
    species%correction_grid=trial%correction_grid
    species%traces=trial%traces
  end subroutine stage_species

  logical function valid_problem_state(problem,state)
    type(actin_sequential_problem_t),intent(in)::problem
    type(actin_sequential_state_t),intent(in)::state
    integer::nx_local,ny_local,nmarker
    valid_problem_state=.false.
    if(.not.allocated(state%network%field) .or. &
         .not.allocated(state%free%field))return
    nx_local=size(state%network%field,1); ny_local=size(state%network%field,2)
    if(any(shape(state%free%field)/=[nx_local,ny_local]))return
    if(.not.allocated(state%network%density) .or. &
         .not.allocated(state%free%density))return
    nmarker=size(state%network%density)
    if(size(state%free%density)/=nmarker .or. nmarker<=0)return
    if(.not.valid_species(state%network,nx_local,ny_local,nmarker))return
    if(.not.valid_species(state%free,nx_local,ny_local,nmarker))return
    if(.not.valid_geometry(state%geometry,nx_local,ny_local,nmarker))return
    if(.not.valid_geometry(problem%next_geometry,nx_local,ny_local,nmarker))return
    if(validate_scalar_operator(problem%network_operator,nx_local,ny_local, &
         nmarker)/=SCALAR_STATUS_OK)return
    if(validate_scalar_operator(problem%free_operator,nx_local,ny_local, &
         nmarker)/=SCALAR_STATUS_OK)return
    if(.not.valid_coefficients(problem%coefficients))return
    if(.not.is_finite_scalar(problem%dt) .or. problem%dt<=0.0_dp)return
    if(.not.allocated(problem%coupling_q))return
    if(size(problem%coupling_q)/=nmarker)return
    if(.not.all(is_finite_scalar(problem%coupling_q)))return
    if(.not.allocated(problem%beta_network) .or. &
         .not.allocated(problem%beta_free))return
    if(.not.allocated(problem%frozen_old_free_trace))return
    if(size(problem%beta_network)/=nmarker .or. &
         size(problem%beta_free)/=nmarker .or. &
         size(problem%frozen_old_free_trace)/=nmarker)return
    if(.not.all(is_finite_scalar(problem%beta_network)) .or. &
         .not.all(is_finite_scalar(problem%beta_free)) .or. &
         .not.all(is_finite_scalar(problem%frozen_old_free_trace)))return
    valid_problem_state=.true.
  end function valid_problem_state

  logical function valid_species(species,nx_local,ny_local,nmarker)
    type(actin_species_state_t),intent(in)::species
    integer,intent(in)::nx_local,ny_local,nmarker
    valid_species=.false.
    if(any(shape(species%field)/=[nx_local,ny_local]))return
    if(.not.all(is_finite_scalar(species%field)))return
    if(.not.all(is_finite_scalar(species%density)))return
    if(.not.allocated(species%correction_coefficients))return
    if(any(shape(species%correction_coefficients)/=[nmarker,6]))return
    if(.not.allocated(species%correction_grid))return
    if(any(shape(species%correction_grid)/=[nx_local,ny_local]))return
    if(.not.all(is_finite_scalar(species%correction_coefficients)))return
    if(.not.all(is_finite_scalar(species%correction_grid)))return
    if(.not.valid_traces(species%traces,nmarker))return
    valid_species=.true.
  end function valid_species

  logical function valid_geometry(geometry,nx_local,ny_local,nmarker)
    type(actin_geometry_history_t),intent(in)::geometry
    integer,intent(in)::nx_local,ny_local,nmarker
    valid_geometry=.false.
    if(.not.allocated(geometry%marker_x) .or. &
         .not.allocated(geometry%marker_y))return
    if(size(geometry%marker_x)/=nmarker .or. &
         size(geometry%marker_y)/=nmarker)return
    if(.not.all(is_finite_scalar(geometry%marker_x)) .or. &
         .not.all(is_finite_scalar(geometry%marker_y)))return
    if(.not.allocated(geometry%nearest_marker))return
    if(any(shape(geometry%nearest_marker)/=[nx_local,ny_local]))return
    if(any(geometry%nearest_marker<1) .or. &
         any(geometry%nearest_marker>nmarker))return
    valid_geometry=.true.
  end function valid_geometry

  logical function valid_trial(trial,nx_local,ny_local,nmarker)
    type(actin_extension_trial_t),intent(in)::trial
    integer,intent(in)::nx_local,ny_local,nmarker
    valid_trial=.false.
    if(.not.allocated(trial%field))return
    if(any(shape(trial%field)/=[nx_local,ny_local]))return
    if(.not.allocated(trial%correction_coefficients))return
    if(any(shape(trial%correction_coefficients)/=[nmarker,6]))return
    if(.not.allocated(trial%correction_grid))return
    if(any(shape(trial%correction_grid)/=[nx_local,ny_local]))return
    if(.not.all(is_finite_scalar(trial%field)))return
    if(.not.all(is_finite_scalar(trial%correction_coefficients)))return
    if(.not.all(is_finite_scalar(trial%correction_grid)))return
    if(.not.valid_traces(trial%traces,nmarker))return
    valid_trial=.true.
  end function valid_trial

  logical function valid_traces(traces,nmarker)
    type(interface_traces_t),intent(in)::traces
    integer,intent(in)::nmarker
    valid_traces=.false.
    if(.not.allocated(traces%value_average))return
    if(size(traces%value_average)/=nmarker)return
    if(.not.allocated(traces%normal_average))return
    if(size(traces%normal_average)/=nmarker)return
    if(.not.allocated(traces%value_interior))return
    if(size(traces%value_interior)/=nmarker)return
    if(.not.allocated(traces%value_exterior))return
    if(size(traces%value_exterior)/=nmarker)return
    if(.not.allocated(traces%normal_interior))return
    if(size(traces%normal_interior)/=nmarker)return
    if(.not.allocated(traces%normal_exterior))return
    if(size(traces%normal_exterior)/=nmarker)return
    if(.not.all(is_finite_scalar(traces%value_average)))return
    if(.not.all(is_finite_scalar(traces%normal_average)))return
    if(.not.all(is_finite_scalar(traces%value_interior)))return
    if(.not.all(is_finite_scalar(traces%value_exterior)))return
    if(.not.all(is_finite_scalar(traces%normal_interior)))return
    if(.not.all(is_finite_scalar(traces%normal_exterior)))return
    valid_traces=.true.
  end function valid_traces

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

  subroutine initialize_result(result)
    type(actin_sequential_result_t),intent(out)::result
    result%status=ACTIN_SEQUENTIAL_INVALID
    result%evaluator_status=0
  end subroutine initialize_result

  logical function final_residual_is_accepted(residual,gmres,options)
    real(dp),intent(in)::residual(:)
    type(gmres_result_t),intent(in)::gmres
    type(gmres_options_t),intent(in)::options
    real(dp)::residual_norm,relative_threshold,threshold
    logical::ok
    final_residual_is_accepted=.false.
    residual_norm=stable_norm(residual,ok); if(.not.ok)return
    call checked_multiply(options%relative_tolerance, &
         gmres%initial_residual_norm,relative_threshold,ok)
    if(.not.ok)relative_threshold=huge(1.0_dp)
    threshold=max(options%absolute_tolerance,relative_threshold, &
         gmres%absolute_residual)+1024.0_dp*epsilon(1.0_dp)
    final_residual_is_accepted=residual_norm<=threshold
  end function final_residual_is_accepted

  real(dp) function stable_norm(values,ok) result(norm_value)
    real(dp),intent(in)::values(:)
    logical,intent(out)::ok
    real(dp)::scale,ssq,magnitude,ratio,root_ssq
    integer::i
    scale=0.0_dp; ssq=1.0_dp; norm_value=0.0_dp; ok=.false.
    if(.not.all(is_finite_scalar(values)))return
    do i=1,size(values)
      magnitude=abs(values(i))
      if(magnitude/=0.0_dp)then
        if(scale<magnitude)then
          ratio=scale/magnitude
          ssq=1.0_dp+ssq*ratio*ratio
          scale=magnitude
        else
          ratio=magnitude/scale
          ssq=ssq+ratio*ratio
        endif
      endif
    enddo
    if(scale==0.0_dp)then; norm_value=0.0_dp; ok=.true.; return; endif
    root_ssq=sqrt(ssq)
    if(scale>huge(1.0_dp)/root_ssq)return
    norm_value=scale*root_ssq
    ok=is_finite_scalar(norm_value)
  end function stable_norm

  subroutine checked_multiply(left,right,value,ok)
    real(dp),intent(in)::left,right
    real(dp),intent(out)::value
    logical,intent(out)::ok
    value=0.0_dp; ok=.false.
    if(.not.is_finite_scalar(left) .or. .not.is_finite_scalar(right))return
    if(left==0.0_dp .or. right==0.0_dp)then; ok=.true.; return; endif
    if(exponent(abs(left))+exponent(abs(right))>=maxexponent(value))return
    value=left*right
    ok=is_finite_scalar(value)
  end subroutine checked_multiply

  subroutine checked_add(left,right,value,ok)
    real(dp),intent(in)::left,right
    real(dp),intent(out)::value
    logical,intent(out)::ok
    value=0.0_dp; ok=.false.
    if(.not.is_finite_scalar(left) .or. .not.is_finite_scalar(right))return
    if(sign(1.0_dp,left)==sign(1.0_dp,right))then
      if(max(exponent(abs(left)),exponent(abs(right)))>= &
           maxexponent(value)-1)return
    endif
    value=left+right
    ok=is_finite_scalar(value)
  end subroutine checked_add

end module fsi_actin_sequential_harness_mod
