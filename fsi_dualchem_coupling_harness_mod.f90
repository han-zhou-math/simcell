! SPDX-License-Identifier: BSD-3-Clause
!
! Minimal Stage-06 one-way coupling harness.
!
! The host FSI solver owns and advances (u,v,p,X).  This module consumes the
! already accepted u^{n+1}, X^n, and X^{n+1}; it never calls the dualchem marker
! advance routine.  Chemical state is evaluated in disposable local copies and
! committed only after BOTH original one-sided solves converge:
!
!     exterior using c_i^n  ->  interior using c_e^{n+1}.
!
! Thus an invalid input or a failed second solve leaves the previously accepted
! chemical fields, interface traces, densities, and correction histories intact.
module fsi_dualchem_coupling_harness_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use parameters, only: dp, nx, ny, npts, xmin, xmax, ymin, ymax, hg, &
      cpi, pi, tupi, scc, dualchem_kc, dualchem_kp, dualchem_diffusion, &
      dualchem_interface_gmres_rtol, &
      dualchem_interface_gmres_atol_physical,dualchem_initial_concentration, &
      dualchem_initial_polarization,zero,one,is_finite_run_scalar
  use fsi_transport_snapshot_mod, only: fsi_transport_snapshot_t, &
      FSI_SNAPSHOT_OK,publish_fsi_transport_snapshot, &
      clear_fsi_transport_snapshot
  use fsi_dualchem_velocity_bridge_mod, only: VELOCITY_BRIDGE_OK, &
      install_dualchem_transport_snapshot,clear_fsi_velocity_bridge, &
      interpolate_fsi_velocity_to_marker,get_fsi_velocity_generation
  use chemical_interface_flux_mod, only: CHEMICAL_EXTERIOR, CHEMICAL_INTERIOR, &
      CHEMICAL_FLUX_INVALID, CHEMICAL_FLUX_OK, chemical_normalization_t, &
      build_chemical_normalization, build_relative_normal_velocity, &
      build_chemical_robin, normalize_chemical_absolute_tolerance
  use grid_types, only: LagrangianGrid, EulerianGrid
  use geometry_mod, only: update_geometry_wrapper
  use linear_solver_mod, only: LinearSolver, jmp_func, &
      evaluate_interface_gmres_residual, GMRES_NOT_CONVERGED
  implicit none
  private

  integer, parameter, public :: DUALCHEM_COUPLING_OK = 0
  integer, parameter, public :: DUALCHEM_COUPLING_INVALID = 1
  integer, parameter, public :: DUALCHEM_COUPLING_NOT_READY = 2
  integer, parameter, public :: DUALCHEM_COUPLING_SOLVE_FAILED = 3

#ifdef SIMCELL_TESTING
  integer, parameter, public :: TASK6_CORRUPT_GMRES_RTOL = 1
  integer, parameter, public :: TASK6_CORRUPT_GMRES_ATOL = 2
  integer, parameter, public :: TASK6_CORRUPT_GMRES_RELATIVE = 3
  integer, save :: task6_gmres_corruption_mode = 0
  integer, save :: task6_gmres_corruption_side = 0
#endif

  type(LagrangianGrid), save :: accepted_lag_in, accepted_lag_out
  type(EulerianGrid), save :: accepted_eul_in, accepted_eul_out
  type(LinearSolver), save :: accepted_solver_in, accepted_solver_out
  real(dp), allocatable, save :: accepted_u_in(:,:), accepted_u_out(:,:)
  real(dp), allocatable, save :: accepted_psi_in(:), accepted_psi_out(:)
  real(dp), allocatable, save :: accepted_trace_in(:), accepted_trace_out(:)
  ! Stage 07 owns the physical difference between the two separately solved
  ! one-sided traces.  It is not lag_grid%phi=[w] (which is zero) and it is
  ! not either derivative-jump density accepted_psi_*.
  real(dp), allocatable, save :: accepted_physical_jump(:)
  ! The accepted snapshot is independent of chemical state and is copied by
  ! value.  Stage 10 will pass this same object to the separate actin harness.
  type(fsi_transport_snapshot_t),save :: accepted_transport_snapshot
  integer,save :: accepted_snapshot_id=0
  logical, save :: coupling_is_initialized = .false.
  integer, save :: accepted_steps = 0

  public :: initialize_dualchem_coupling
  public :: advance_dualchem_oneway
  public :: get_dualchem_mass
  public :: get_accepted_concentration_jump
  public :: get_accepted_physical_chemical_field
  public :: get_accepted_marker_geometry
  public :: get_accepted_transport_snapshot
  public :: reconstruct_physical_concentration_jump
  public :: finalize_dualchem_coupling
#ifdef SIMCELL_TESTING
  public :: arm_task6_gmres_record_corruption
#endif

contains

#ifdef SIMCELL_TESTING
  subroutine arm_task6_gmres_record_corruption(mode,side,status)
    integer,intent(in) :: mode,side
    integer,intent(out) :: status

    status=DUALCHEM_COUPLING_INVALID
    if(mode<TASK6_CORRUPT_GMRES_RTOL .or. &
       mode>TASK6_CORRUPT_GMRES_RELATIVE)return
    if(side/=CHEMICAL_EXTERIOR .and. side/=CHEMICAL_INTERIOR)return
    task6_gmres_corruption_mode=mode
    task6_gmres_corruption_side=side
    status=DUALCHEM_COUPLING_OK
  end subroutine arm_task6_gmres_record_corruption

  subroutine apply_task6_gmres_record_corruption(solver,side)
    type(LinearSolver),intent(inout) :: solver
    integer,intent(in) :: side

    if(task6_gmres_corruption_side/=side)return
    select case(task6_gmres_corruption_mode)
    case(TASK6_CORRUPT_GMRES_RTOL)
      solver%tol=2.0_dp*dualchem_interface_gmres_rtol
    case(TASK6_CORRUPT_GMRES_ATOL)
      solver%gmres_atol_normalized=0.5_dp*solver%gmres_atol_normalized
    case(TASK6_CORRUPT_GMRES_RELATIVE)
      solver%gmres_final_relative_residual=zero
    case default
      return
    end select
    task6_gmres_corruption_mode=0
    task6_gmres_corruption_side=0
  end subroutine apply_task6_gmres_record_corruption
#endif

  subroutine initialize_dualchem_coupling(x,y,status)
    real(dp), intent(in) :: x(:),y(:)
    integer, intent(out) :: status
    integer :: i,j
    real(dp) :: cell_x,polarized_value,internal_baseline

    call finalize_dualchem_coupling()
    status=DUALCHEM_COUPLING_INVALID
    if(size(x)/=npts .or. size(y)/=npts)return
    if(.not.all(ieee_is_finite(x)) .or. .not.all(ieee_is_finite(y)))return
    if(.not.markers_are_wall_safe(x,y))return
    ! Direct callers must satisfy the same run-level contract as fmain/readpar.
    ! Reject before chemical grids, solvers, or state arrays are initialized.
    if(.not.is_finite_run_scalar(dualchem_diffusion))return
    if(dualchem_diffusion<=zero)return

    ! Starter normally initializes pi/tupi in InitLocCell.  Setting the shared
    ! values here as well keeps this harness independently testable without
    ! creating a second parameter module or a hidden dualchem parameter copy.
    pi=cpi
    tupi=2.0_dp*cpi
    do i=1,npts
      scc(i)=real(i-1,dp)*tupi/real(npts,dp)
    end do

    call accepted_lag_in%init(npts)
    call accepted_lag_out%init(npts)
    call accepted_eul_in%init(nx,ny,xmin,xmax,ymin,ymax)
    call accepted_eul_out%init(nx,ny,xmin,xmax,ymin,ymax)
    accepted_lag_in%x=x
    accepted_lag_in%y=y
    accepted_lag_in%x_old=x
    accepted_lag_in%y_old=y
    call update_geometry_wrapper(accepted_lag_in,accepted_eul_in,one,-1)
    call copy_lagrangian_geometry(accepted_lag_in,accepted_lag_out)
    call copy_eulerian_geometry(accepted_eul_in,accepted_eul_out)

    call accepted_solver_in%init(100,50,dualchem_interface_gmres_rtol,npts)
    call accepted_solver_out%init(100,50,dualchem_interface_gmres_rtol,npts)
    accepted_solver_in%bc_n=zero
    accepted_solver_in%bc_s=zero
    accepted_solver_out%bc_n=zero
    accepted_solver_out%bc_s=zero

    allocate(accepted_u_in(nx,ny),accepted_u_out(nx,ny))
    allocate(accepted_psi_in(npts),accepted_psi_out(npts))
    allocate(accepted_trace_in(npts),accepted_trace_out(npts))
    allocate(accepted_physical_jump(npts))
    ! The transport solver advances w=c-1.  Keep runtime input in terms of the
    ! physical dimensionless concentration c, so c=1 initializes w=0.
    internal_baseline=dualchem_initial_concentration-one
    do j=1,ny
      do i=1,nx
        cell_x=xmin+(real(i,dp)-0.5_dp)*hg
        polarized_value=dualchem_initial_polarization*cell_x
        accepted_u_in(i,j)=internal_baseline+polarized_value
        accepted_u_out(i,j)=internal_baseline-polarized_value
      enddo
    enddo
    accepted_psi_in=zero
    accepted_psi_out=zero
    accepted_trace_in=internal_baseline+dualchem_initial_polarization*x
    accepted_trace_out=internal_baseline-dualchem_initial_polarization*x
    accepted_physical_jump=accepted_trace_in-accepted_trace_out
    accepted_steps=0
    coupling_is_initialized=.true.
    status=DUALCHEM_COUPLING_OK
  end subroutine initialize_dualchem_coupling

  subroutine advance_dualchem_oneway(u,v,x_old,y_old,x_new,y_new,dt,time,status)
    real(dp), intent(in) :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1)
    real(dp), intent(in) :: x_old(:),y_old(:),x_new(:),y_new(:)
    real(dp), intent(in) :: dt,time
    integer, intent(out) :: status

    type(LagrangianGrid) :: lag_in,lag_out
    type(EulerianGrid) :: eul_in,eul_out
    type(LinearSolver) :: solver_in,solver_out
    type(fsi_transport_snapshot_t),allocatable :: candidate_snapshot
    real(dp) :: psi_in(npts),psi_out(npts),trace_in(npts),trace_out(npts)
    real(dp) :: candidate_physical_jump(npts)
    real(dp) :: fjmp(npts),zeros(npts),h0(npts),rhs(npts),work(npts)
    real(dp) :: marker_velocity(npts,2),v_gamma(npts,2),r(npts)
    real(dp) :: v_gamma_check(npts,2),r_check(npts)
    real(dp) :: unknown_coef,gradient_coef,opposite_coef,shift_constant
    real(dp) :: abs_resid_in,rel_resid_in,abs_resid_out,rel_resid_out
    real(dp) :: physical_abs_resid_in,physical_abs_resid_out
    logical :: first_step,converged,gate_accepted
    integer :: i,isel,bridge_status,snapshot_status,allocation_status
    integer :: flux_status
    integer :: candidate_snapshot_id
    integer(int64) :: generation,generation_check
    type(chemical_normalization_t) :: normalization

    status=DUALCHEM_COUPLING_INVALID
    if(.not.coupling_is_initialized)then
      status=DUALCHEM_COUPLING_NOT_READY
      return
    end if
    if(size(x_old)/=npts .or. size(y_old)/=npts .or. &
       size(x_new)/=npts .or. size(y_new)/=npts)return
    if(.not.ieee_is_finite(dt) .or. dt<=zero)return
    if(.not.ieee_is_finite(time))return
    ! Accept and freeze the one common D before any trial geometry, snapshot,
    ! chemical allocation, or affine-operator evaluation can occur.
    call build_chemical_normalization(dualchem_diffusion,dt,normalization, &
        flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK)return
    if(.not.all(ieee_is_finite(x_old)) .or. .not.all(ieee_is_finite(y_old)))return
    if(.not.all(ieee_is_finite(x_new)) .or. .not.all(ieee_is_finite(y_new)))return
    if(.not.markers_are_wall_safe(x_old,y_old))return
    if(.not.markers_are_wall_safe(x_new,y_new))return

    ! Intrinsic assignment deep-copies every allocatable component.  These
    ! local objects are the trial transaction; accepted module state below is
    ! untouched until the final commit block.
    lag_in=accepted_lag_in
    lag_out=accepted_lag_out
    eul_in=accepted_eul_in
    eul_out=accepted_eul_out
    solver_in=accepted_solver_in
    solver_out=accepted_solver_out
    psi_in=accepted_psi_in
    psi_out=accepted_psi_out
    trace_in=accepted_trace_in
    trace_out=accepted_trace_out

    lag_in%x_old=x_old
    lag_in%y_old=y_old
    lag_in%x=x_new
    lag_in%y=y_new
    first_step=(accepted_steps==0)
    isel=merge(-1,0,first_step)
    call update_geometry_wrapper(lag_in,eul_in,dt,isel)
    call copy_lagrangian_geometry(lag_in,lag_out)
    call copy_eulerian_geometry(eul_in,eul_out)

    ! Publish exactly one complete accepted-FSI candidate after the current
    ! geometry has supplied its outward normals.  The dualchem adapter installs
    ! this value object before H(0); every GMRES callback therefore sees the
    ! same face and marker velocity, geometry history, dt, and time level.
    allocate(candidate_snapshot,stat=allocation_status)
    if(allocation_status/=0)return
    candidate_snapshot_id=accepted_steps+1
    call publish_fsi_transport_snapshot(candidate_snapshot,u,v,x_old,y_old, &
        x_new,y_new,lag_in%normal,dt,time,candidate_snapshot_id, &
        snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)return
    call install_dualchem_transport_snapshot(candidate_snapshot, &
        candidate_snapshot_id,bridge_status)
    if(bridge_status/=VELOCITY_BRIDGE_OK)then
      call restore_accepted_transport_snapshot()
      return
    end if

    ! Stage 6A moving-interface flux correction.  The installed snapshot already
    ! holds the accepted u^{n+1} faces; interpolate the marker velocity at
    ! X^{n+1}, form the interface velocity V_Gamma=(X^{n+1}-X^n)/dt, and derive
    ! the relative normal velocity r=(v_marker-V_Gamma).n once for both
    ! one-sided solves.  The installed snapshot id is the bridge generation the
    ! affine operator will later require to match.
    do i=1,npts
      call interpolate_fsi_velocity_to_marker(x_new(i),y_new(i), &
          marker_velocity(i,1),marker_velocity(i,2),bridge_status)
      if(bridge_status/=VELOCITY_BRIDGE_OK)then
        call restore_accepted_transport_snapshot()
        return
      end if
      v_gamma(i,:)=[(x_new(i)-x_old(i))/dt,(y_new(i)-y_old(i))/dt]
    end do
    call build_relative_normal_velocity(marker_velocity,v_gamma,lag_in%normal, &
        r,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK)then
      call restore_accepted_transport_snapshot()
      return
    end if
    call get_fsi_velocity_generation(generation,bridge_status)
    if(bridge_status/=VELOCITY_BRIDGE_OK)then
      call restore_accepted_transport_snapshot()
      return
    end if
    call solver_out%set_interface_kinematics(v_gamma,r,generation,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK)then
      call restore_accepted_transport_snapshot()
      return
    end if
    call solver_in%set_interface_kinematics(v_gamma,r,generation,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK)then
      call restore_accepted_transport_snapshot()
      return
    end if
    ! Round-trip gate: the frozen kinematics and the bridge generation must
    ! agree on both solvers before any affine-operator evaluation may consume
    ! them.  A mismatch here is a stale trial, never a usable Krylov vector.
    call get_fsi_velocity_generation(generation_check,bridge_status)
    if(bridge_status/=VELOCITY_BRIDGE_OK .or. generation_check/=generation)then
      call restore_accepted_transport_snapshot()
      return
    end if
    call solver_out%get_interface_kinematics(v_gamma_check,r_check, &
        generation_check,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK .or. generation_check/=generation .or. &
       maxval(abs(v_gamma_check-v_gamma))>zero .or. &
       maxval(abs(r_check-r))>zero)then
      call restore_accepted_transport_snapshot()
      return
    end if
    call solver_in%get_interface_kinematics(v_gamma_check,r_check, &
        generation_check,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK .or. generation_check/=generation .or. &
       maxval(abs(v_gamma_check-v_gamma))>zero .or. &
       maxval(abs(r_check-r))>zero)then
      call restore_accepted_transport_snapshot()
      return
    end if

    ! Exact invariant state used by the JTB actin-only configurations.  With
    ! zero membrane permeability, zero pump, and no initial polarization, the
    ! two chemical fields are the same spatial constant.  Advection and
    ! diffusion preserve that constant exactly, so an auxiliary interface
    ! GMRES solve has a numerically zero right-hand side and can only measure
    ! matrix-free roundoff.  Commit the updated geometry/snapshot directly;
    ! active-chemistry configurations continue through the normal solve below.
    if(dualchem_kc==zero .and. dualchem_kp==zero .and. &
         dualchem_initial_polarization==zero)then
      accepted_lag_in=lag_in
      accepted_lag_out=lag_out
      accepted_eul_in=eul_in
      accepted_eul_out=eul_out
      accepted_solver_in=solver_in
      accepted_solver_out=solver_out
      accepted_physical_jump=zero
      accepted_transport_snapshot=candidate_snapshot
      accepted_snapshot_id=candidate_snapshot_id
      accepted_steps=accepted_steps+1
      status=DUALCHEM_COUPLING_OK
      call report_chemical_solver(candidate_snapshot_id,CHEMICAL_EXTERIOR, &
           1,solver_out,zero)
      call report_chemical_solver(candidate_snapshot_id,CHEMICAL_INTERIOR, &
           1,solver_in,zero)
      return
    end if

    zeros=zero
    fjmp=zero

    ! Exterior auxiliary solve, coupled to the previously accepted interior
    ! trace.  Stored fields are shifted by c0=1, exactly as the original driver.
    call solver_out%set_u_prev(accepted_u_out,dt)
    solver_out%bc_n=zero
    solver_out%bc_s=zero
    solver_out%eulerian_field=zero
    lag_out%phi=zero
    do i=1,npts
      call build_chemical_robin(CHEMICAL_EXTERIOR,r(i),dualchem_kc, &
          jmp_func(scc(i),time),normalization%diffusion,unknown_coef, &
          gradient_coef,opposite_coef,shift_constant,flux_status)
      if(flux_status/=CHEMICAL_FLUX_OK)then
        call restore_accepted_transport_snapshot()
        return
      end if
      lag_out%rhs(i)=opposite_coef*accepted_trace_in(i)+shift_constant
    end do
    call solver_out%compute_H(lag_out,eul_out,zeros,h0,zero,dt,time,isel, &
        first_step,-1,fjmp,.false.,.false.)
    if(solver_out%inner_solve_failed)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if
    rhs=-h0
    call solver_out%solve(lag_out,eul_out,rhs,psi_out,zero,dt,time,isel, &
        first_step,-1,fjmp,.false.,.false.,converged_out=converged)
#ifdef SIMCELL_TESTING
    call apply_task6_gmres_record_corruption(solver_out,CHEMICAL_EXTERIOR)
#endif
    call independently_validate_gmres_record(solver_out,converged, &
        normalization%diffusion,gate_accepted,abs_resid_out,rel_resid_out)
    if(.not.gate_accepted)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if
    call solver_out%compute_H(lag_out,eul_out,psi_out,work,zero,dt,time,isel, &
        first_step,-1,fjmp,.false.,.false.,u_at_interface_out=trace_out)
    if(solver_out%inner_solve_failed)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if

    ! Interior auxiliary solve uses the newly converged exterior trace, which
    ! preserves the original sequential (Gauss-Seidel-like) split.
    call solver_in%set_u_prev(accepted_u_in,dt)
    solver_in%bc_n=zero
    solver_in%bc_s=zero
    solver_in%eulerian_field=zero
    lag_in%phi=zero
    do i=1,npts
      call build_chemical_robin(CHEMICAL_INTERIOR,r(i),dualchem_kc, &
          jmp_func(scc(i),time),normalization%diffusion,unknown_coef, &
          gradient_coef,opposite_coef,shift_constant,flux_status)
      if(flux_status/=CHEMICAL_FLUX_OK)then
        call restore_accepted_transport_snapshot()
        return
      end if
      lag_in%rhs(i)=opposite_coef*trace_out(i)+shift_constant
    end do
    call solver_in%compute_H(lag_in,eul_in,zeros,h0,zero,dt,time,isel, &
        first_step,1,fjmp,.false.,.false.)
    if(solver_in%inner_solve_failed)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if
    rhs=-h0
    call solver_in%solve(lag_in,eul_in,rhs,psi_in,zero,dt,time,isel, &
        first_step,1,fjmp,.false.,.false.,converged_out=converged)
#ifdef SIMCELL_TESTING
    call apply_task6_gmres_record_corruption(solver_in,CHEMICAL_INTERIOR)
#endif
    call independently_validate_gmres_record(solver_in,converged, &
        normalization%diffusion,gate_accepted,abs_resid_in,rel_resid_in)
    if(.not.gate_accepted)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if
    call solver_in%compute_H(lag_in,eul_in,psi_in,work,zero,dt,time,isel, &
        first_step,1,fjmp,.false.,.false.,u_at_interface_in=trace_in)
    if(solver_in%inner_solve_failed)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if

    ! Stage 6A transaction gate: the trial is acceptable only when the
    ! recomputed (not recursive) GMRES residual of BOTH one-sided solves meets
    ! the declared interface tolerance.
    call convert_absolute_residual_to_physical(normalization%diffusion, &
        abs_resid_in,physical_abs_resid_in,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK .or. physical_abs_resid_in<zero)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if
    call convert_absolute_residual_to_physical(normalization%diffusion, &
        abs_resid_out,physical_abs_resid_out,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK .or. physical_abs_resid_out<zero)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if

    if(.not.all(ieee_is_finite(eul_in%u(1:nx,1:ny))) .or. &
       .not.all(ieee_is_finite(eul_out%u(1:nx,1:ny))) .or. &
       .not.all(ieee_is_finite(trace_in)) .or. &
       .not.all(ieee_is_finite(trace_out)))then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if

    ! The physical membrane feedback is reconstructed from the two physical
    ! one-sided values.  Both fields are stored as w=c-1, so the common shift
    ! cancels exactly: (w_i+1)-(w_e+1)=w_i-w_e.  This quantity is unrelated
    ! to the auxiliary correction jump [w]=0 used inside each one-sided solve.
    call reconstruct_physical_concentration_jump(trace_in,trace_out, &
        candidate_physical_jump,status)
    if(status/=DUALCHEM_COUPLING_OK)then
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      call restore_accepted_transport_snapshot()
      return
    end if

    ! Atomic commit: no earlier branch writes accepted chemical state.
    call solver_in%update_history()
    call solver_out%update_history()
    accepted_lag_in=lag_in
    accepted_lag_out=lag_out
    accepted_eul_in=eul_in
    accepted_eul_out=eul_out
    accepted_solver_in=solver_in
    accepted_solver_out=solver_out
    accepted_u_in=eul_in%u(1:nx,1:ny)
    accepted_u_out=eul_out%u(1:nx,1:ny)
    accepted_psi_in=psi_in
    accepted_psi_out=psi_out
    accepted_trace_in=trace_in
    accepted_trace_out=trace_out
    accepted_physical_jump=candidate_physical_jump
    accepted_transport_snapshot=candidate_snapshot
    accepted_snapshot_id=candidate_snapshot_id
    accepted_steps=accepted_steps+1
    status=DUALCHEM_COUPLING_OK
    call report_chemical_solver(candidate_snapshot_id,CHEMICAL_EXTERIOR,0, &
         solver_out,physical_abs_resid_out)
    call report_chemical_solver(candidate_snapshot_id,CHEMICAL_INTERIOR,0, &
         solver_in,physical_abs_resid_in)
  end subroutine advance_dualchem_oneway

  subroutine report_chemical_solver(snapshot_id,side,bypass,solver, &
      physical_absolute_residual)
    integer,intent(in)::snapshot_id,side,bypass
    type(LinearSolver),intent(in)::solver
    real(dp),intent(in)::physical_absolute_residual
    integer::converged,inner_failed,reason,iterations
    real(dp)::initial_residual,normalized_absolute_residual
    real(dp)::relative_residual,rtol,physical_atol

    converged=0; inner_failed=0; reason=0; iterations=0
    initial_residual=zero; normalized_absolute_residual=zero
    relative_residual=zero; rtol=dualchem_interface_gmres_rtol
    physical_atol=dualchem_interface_gmres_atol_physical
    if(bypass==0)then
      converged=merge(1,0,solver%gmres_converged)
      inner_failed=merge(1,0,solver%inner_solve_failed)
      reason=solver%gmres_convergence_reason
      iterations=solver%gmres_total_krylov_iterations
      initial_residual=solver%gmres_initial_normalized_residual
      normalized_absolute_residual=solver%gmres_final_normalized_residual
      relative_residual=solver%gmres_final_relative_residual
      rtol=solver%tol
    endif
    write(*,'(a,1x,8(i0,1x),6(es24.16,1x))') &
         'FIG2_CHEMICAL_SOLVER',snapshot_id,side,bypass, &
         DUALCHEM_COUPLING_OK,converged,reason,iterations,inner_failed, &
         initial_residual,normalized_absolute_residual, &
         physical_absolute_residual,relative_residual,rtol,physical_atol
  end subroutine report_chemical_solver

  subroutine get_dualchem_mass(total_mass,status)
    real(dp), intent(out) :: total_mass
    integer, intent(out) :: status
    integer :: i,j

    total_mass=zero
    if(.not.coupling_is_initialized)then
      status=DUALCHEM_COUPLING_NOT_READY
      return
    end if
    do j=1,ny
      do i=1,nx
        if(accepted_eul_in%idf(i,j)==1)then
          total_mass=total_mass+(accepted_u_in(i,j)+one)*hg*hg
        else if(accepted_eul_out%idf(i,j)==-1)then
          total_mass=total_mass+(accepted_u_out(i,j)+one)*hg*hg
        end if
      end do
    end do
    if(.not.ieee_is_finite(total_mass))then
      total_mass=zero
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      return
    end if
    status=DUALCHEM_COUPLING_OK
  end subroutine get_dualchem_mass

  subroutine reconstruct_physical_concentration_jump(interior_trace, &
      exterior_trace,jump,status)
    real(dp), intent(in) :: interior_trace(:),exterior_trace(:)
    real(dp), intent(out) :: jump(:)
    integer, intent(out) :: status

    jump=zero
    status=DUALCHEM_COUPLING_INVALID
    if(size(interior_trace)/=npts .or. size(exterior_trace)/=npts .or. &
       size(jump)/=npts)return
    if(.not.all(ieee_is_finite(interior_trace)))return
    if(.not.all(ieee_is_finite(exterior_trace)))return
    jump=interior_trace-exterior_trace
    if(.not.all(ieee_is_finite(jump)))then
      jump=zero
      return
    end if
    status=DUALCHEM_COUPLING_OK
  end subroutine reconstruct_physical_concentration_jump

  subroutine get_accepted_concentration_jump(jump,status)
    real(dp), intent(out) :: jump(:)
    integer, intent(out) :: status

    jump=zero
    if(.not.coupling_is_initialized)then
      status=DUALCHEM_COUPLING_NOT_READY
      return
    end if
    if(size(jump)/=npts)then
      status=DUALCHEM_COUPLING_INVALID
      return
    end if
    jump=accepted_physical_jump
    status=DUALCHEM_COUPLING_OK
  end subroutine get_accepted_concentration_jump

  subroutine get_accepted_physical_chemical_field(field,status)
    real(dp), intent(out) :: field(:,:)
    integer, intent(out) :: status
    integer :: i,j

    ! Refinement/audit view of the physical piecewise field.  The auxiliary
    ! full-box extensions remain private: at each cell center select the
    ! accepted solution on the physical side and undo the common w=c-1 shift.
    field=zero
    if(.not.coupling_is_initialized)then
      status=DUALCHEM_COUPLING_NOT_READY
      return
    end if
    if(size(field,1)/=nx .or. size(field,2)/=ny)then
      status=DUALCHEM_COUPLING_INVALID
      return
    end if
    do j=1,ny
      do i=1,nx
        if(accepted_eul_in%idf(i,j)==1)then
          field(i,j)=accepted_u_in(i,j)+one
        else if(accepted_eul_out%idf(i,j)==-1)then
          field(i,j)=accepted_u_out(i,j)+one
        else
          field=zero
          status=DUALCHEM_COUPLING_SOLVE_FAILED
          return
        end if
      end do
    end do
    if(.not.all(ieee_is_finite(field)))then
      field=zero
      status=DUALCHEM_COUPLING_SOLVE_FAILED
      return
    end if
    status=DUALCHEM_COUPLING_OK
  end subroutine get_accepted_physical_chemical_field

  subroutine get_accepted_marker_geometry(x,y,status)
    real(dp), intent(out) :: x(:),y(:)
    integer, intent(out) :: status

    ! Read-only audit seam.  The accepted chemical geometry is private so a
    ! caller cannot accidentally alter the transactional state.  Stage 07 uses
    ! this getter to verify that a rejected chemical step changes neither the
    ! accepted concentration jump nor the geometry on which it is defined.
    x=zero
    y=zero
    if(.not.coupling_is_initialized)then
      status=DUALCHEM_COUPLING_NOT_READY
      return
    end if
    if(size(x)/=npts .or. size(y)/=npts)then
      status=DUALCHEM_COUPLING_INVALID
      return
    end if
    x=accepted_lag_in%x
    y=accepted_lag_in%y
    status=DUALCHEM_COUPLING_OK
  end subroutine get_accepted_marker_geometry

  subroutine get_accepted_transport_snapshot(snapshot,snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(out)::snapshot
    integer,intent(out)::snapshot_id,status

    call clear_fsi_transport_snapshot(snapshot)
    snapshot_id=0
    if(.not.coupling_is_initialized .or. accepted_snapshot_id<=0)then
      status=DUALCHEM_COUPLING_NOT_READY
      return
    end if
    snapshot=accepted_transport_snapshot
    snapshot_id=accepted_snapshot_id
    status=DUALCHEM_COUPLING_OK
  end subroutine get_accepted_transport_snapshot

  subroutine restore_accepted_transport_snapshot()
    integer::bridge_status
    if(accepted_snapshot_id<=0)then
      call clear_fsi_velocity_bridge()
      return
    end if
    call install_dualchem_transport_snapshot(accepted_transport_snapshot, &
        accepted_snapshot_id,bridge_status)
    if(bridge_status/=VELOCITY_BRIDGE_OK)call clear_fsi_velocity_bridge()
  end subroutine restore_accepted_transport_snapshot

  subroutine independently_validate_gmres_record(solver,solve_converged, &
      validated_diffusion,accepted,stored_absolute,stored_relative)
    type(LinearSolver),intent(in) :: solver
    logical,intent(in) :: solve_converged
    real(dp),intent(in) :: validated_diffusion
    logical,intent(out) :: accepted
    real(dp),intent(out) :: stored_absolute,stored_relative
    real(dp) :: expected_atol_normalized,gate_relative
    logical :: gate_converged,record_matches
    integer :: residual_status,atol_status,gate_reason

    accepted=.false.
    stored_absolute=huge(one)
    stored_relative=huge(one)
    call normalize_chemical_absolute_tolerance( &
        dualchem_interface_gmres_atol_physical,validated_diffusion, &
        expected_atol_normalized,atol_status)
    if(atol_status/=CHEMICAL_FLUX_OK)return
    call solver%get_last_true_residual(stored_absolute,stored_relative, &
        residual_status)
    call evaluate_interface_gmres_residual(stored_absolute, &
        solver%gmres_initial_normalized_residual, &
        dualchem_interface_gmres_rtol,expected_atol_normalized, &
        gate_converged,gate_reason,gate_relative)

    record_matches= &
        backend_reals_bit_equal(solver%tol,dualchem_interface_gmres_rtol) .and. &
        backend_reals_bit_equal(solver%gmres_atol_normalized, &
        expected_atol_normalized) .and. &
        backend_reals_bit_equal(stored_relative,gate_relative) .and. &
        backend_reals_bit_equal(solver%gmres_diffusion,validated_diffusion)
    accepted=solve_converged .and. residual_status==CHEMICAL_FLUX_OK .and. &
        gate_converged .and. gate_reason/=GMRES_NOT_CONVERGED .and. &
        record_matches .and. &
        gate_reason==solver%gmres_convergence_reason .and. &
        (gate_converged.eqv.solver%gmres_converged)
  end subroutine independently_validate_gmres_record

  pure logical function backend_reals_bit_equal(left,right) result(equal)
    real(dp),intent(in) :: left,right
    equal=transfer(left,0_int64)==transfer(right,0_int64)
  end function backend_reals_bit_equal

  subroutine finalize_dualchem_coupling()
    if(coupling_is_initialized)then
      call accepted_solver_in%clean()
      call accepted_solver_out%clean()
      call accepted_lag_in%clean()
      call accepted_lag_out%clean()
      call accepted_eul_in%clean()
      call accepted_eul_out%clean()
    end if
    if(allocated(accepted_u_in))deallocate(accepted_u_in)
    if(allocated(accepted_u_out))deallocate(accepted_u_out)
    if(allocated(accepted_psi_in))deallocate(accepted_psi_in)
    if(allocated(accepted_psi_out))deallocate(accepted_psi_out)
    if(allocated(accepted_trace_in))deallocate(accepted_trace_in)
    if(allocated(accepted_trace_out))deallocate(accepted_trace_out)
    if(allocated(accepted_physical_jump))deallocate(accepted_physical_jump)
    call clear_fsi_velocity_bridge()
    call clear_fsi_transport_snapshot(accepted_transport_snapshot)
    accepted_snapshot_id=0
    accepted_steps=0
    coupling_is_initialized=.false.
#ifdef SIMCELL_TESTING
    task6_gmres_corruption_mode=0
    task6_gmres_corruption_side=0
#endif
  end subroutine finalize_dualchem_coupling

  logical function markers_are_wall_safe(x,y)
    real(dp), intent(in) :: x(:),y(:)
    markers_are_wall_safe=all(x>=xmin .and. x<=xmax) .and. &
        all(y>=ymin+2.0_dp*hg .and. y<=ymax-2.0_dp*hg)
  end function markers_are_wall_safe

  pure subroutine convert_absolute_residual_to_physical(diffusion,normalized, &
      physical,status)
    real(dp),intent(in) :: diffusion,normalized
    real(dp),intent(out) :: physical
    integer,intent(out) :: status

    physical=zero
    status=CHEMICAL_FLUX_INVALID
    if(.not.ieee_is_finite(diffusion) .or. diffusion<=zero)return
    if(.not.ieee_is_finite(normalized) .or. normalized<zero)return
    if(diffusion>one)then
      if(normalized>huge(one)/diffusion)return
    endif
    physical=diffusion*normalized
    if(.not.ieee_is_finite(physical) .or. physical<zero)then
      physical=zero
      return
    endif
    status=CHEMICAL_FLUX_OK
  end subroutine convert_absolute_residual_to_physical

  subroutine copy_lagrangian_geometry(source,destination)
    type(LagrangianGrid), intent(in) :: source
    type(LagrangianGrid), intent(inout) :: destination
    destination%x=source%x
    destination%y=source%y
    destination%x_old=source%x_old
    destination%y_old=source%y_old
    destination%normal=source%normal
    destination%tangent=source%tangent
    destination%mk=source%mk
    destination%plinkij=source%plinkij
  end subroutine copy_lagrangian_geometry

  subroutine copy_eulerian_geometry(source,destination)
    type(EulerianGrid), intent(in) :: source
    type(EulerianGrid), intent(inout) :: destination
    destination%id=source%id
    destination%idf=source%idf
    destination%idn=source%idn
    destination%dmapc=source%dmapc
    destination%dmapo=source%dmapo
    destination%kaic=source%kaic
    destination%kaio=source%kaio
    destination%chkero=source%chkero
    destination%chkerc=source%chkerc
    destination%oid=source%oid
    destination%qid=source%qid
  end subroutine copy_eulerian_geometry

end module fsi_dualchem_coupling_harness_mod
