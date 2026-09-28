!---------------------------------------------------------------------------  
!Main program of standard IBM
! MODULE: main program 
!> @author 
!> Lingxing Yao
!
!DESCRIPTION: 
!@brief
!> main program for solving Stokes flow with IBM, fixed wall at top/bottom
!> and left/right periodic
!> 
!> NOTE: FORTRAN DOES NOT DISTINGUISH UPPER AND LOWER CASE OF VARIABLES
!--------------------------------------------------------------------------- 
!
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#ifndef PetscCallA
#define PetscCallA(a) call a
#endif
#ifndef PetscCallMPIA
#define PetscCallMPIA(a) call a
#endif
#include <petsc/finclude/petscksp.h>
!
!module linmatop
!#include "petsc/finclude/petscksp.h"
!#include "petsc/finclude/petscvec.h"
!  use petscksp
!  use petscmat
!  use petscvec
!  use petscsys
!
!contains
!end module linmatop
!
program IBM
  use small_solver_mod, only: report_small_solver_audit
  use geometry_mod, only: report_fixed_stencil_audit
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!  use petsc
  use petscsys
  use petscvec

  use parameters
  use geometry
  use linsys
  use myfft
  use IBmod
  use IBforce
  use FSIsolve
  ! Stage 07 keeps the Stage-06 post-FSI chemical transaction and adds exactly
  ! one reverse datum: the previously accepted physical jump c_i-c_e.  No
  ! current-step chemical iterate enters fsisolve.
  use fsi_dualchem_coupling_harness_mod, only: DUALCHEM_COUPLING_OK, &
      initialize_dualchem_coupling, advance_dualchem_oneway, &
      get_dualchem_mass, get_accepted_concentration_jump, &
      get_accepted_physical_chemical_field, &
      get_accepted_transport_snapshot, &
      finalize_dualchem_coupling
  use fsi_transport_snapshot_mod, only: fsi_transport_snapshot_t
  use actin_model_mod, only: actin_coefficients_t,ACTIN_OK, &
      initialize_actin_coefficients
  use fsi_actin_oneway_harness_mod, only: ACTIN_ONEWAY_OK, &
      actin_oneway_manager_t,initialize_actin_oneway,advance_actin_oneway, &
      get_actin_oneway_mass,get_actin_oneway_fields, &
      get_actin_oneway_feedback_state
  use fsi_actin_feedback_mod, only: ACTIN_FEEDBACK_OK, &
      actin_fsi_feedback_t,build_actin_fsi_feedback
!  use chemical_mod
!  use network_mod
!  use solver_mod
!  use energy
!  use linmatop

  implicit none
!
!========================================================================
  double precision :: t0, tf, dt
  double precision :: Esol, Emem, Iflow, En, Inet, Imem, Jmem, Jmem1
!
  integer :: it, ik,jk, Nout, kin, kout, snapshot_unit, snapshot_iostat
  integer :: substep_count,split_count
!
  double precision, dimension(-1:nx+1,-1:ny+1) :: mb1, ma1, &
    u1, u0, v1, v0, ua, va, pa, p0, p1, tiu, tiv, frs, grs
  double precision, dimension(-1:nx+1,-1:ny+1) :: accepted_u,accepted_v, &
    accepted_p
!rm  double precision, dimension(-1:nx+1,-1:ny+1) :: thetu,thetv, sigx,sigy,&
!rm          thetx,thety, grsgx,grsgy,grx, gry, sga
  double precision, dimension(mpts) :: fsx,fsy, xpa,ypa, jv, vbx,vby, usx,usy, vsx,vsy, nvx, nvy
  double precision,dimension(mpts)::accepted_xpt,accepted_ypt
  double precision :: time, time1, time2, time3!, gp(mpts), us(mpts,2)!, tus(mpts,2)
  integer :: ins!, i, N, M
  double precision ::  tp1, tp2, tp3, tp4, tmp, ds
!========================================================================
  double precision :: alphab, alpha, beta, gam
  double precision, dimension(-1:nx+1,4) :: uvbc0!, uvbc1, uvbca !ubt, utp, vbt, vtp
!========================================================================
  integer :: iflg, strlen
  integer :: dualchem_status
  double precision :: dualchem_mass
  double precision, dimension(mpts) :: chemical_jump
  double precision, dimension(nx,ny) :: physical_chemical_field
  type(fsi_transport_snapshot_t) :: accepted_transport_snapshot
  type(actin_coefficients_t) :: actin_coefficients
  type(actin_oneway_manager_t) :: actin_manager
  type(actin_fsi_feedback_t)::actin_feedback
  integer :: actin_status,actin_feedback_status,actin_source_id, &
      fsi_actin_source_id,shared_snapshot_id
  integer :: brinkman_generation,brinkman_solve_calls, &
      brinkman_inner_iterations
  integer :: fsi_solver_status,fsi_outer_iterations, &
      brinkman_krylov_calls,brinkman_worst_reason
  double precision :: fsi_initial_residual,fsi_true_absolute_residual, &
      fsi_true_relative_residual,brinkman_max_true_absolute_residual, &
      brinkman_max_true_relative_residual
  double precision :: actin_mass,max_actin_drag,max_actin_bulk_force, &
      max_actin_stress,step_x_shift,cumulative_x_shift,accepted_time, &
      remaining_time,trial_dt,max_normal_displacement,max_normal_velocity, &
      normal_cfl
  integer,parameter::max_normal_step_splits=50
  double precision, dimension(nx,ny) :: physical_network_actin, &
      physical_free_actin,accepted_network_auxiliary, &
      accepted_network_correction
  integer,dimension(nx,ny)::accepted_network_cell_side
  double precision,dimension(mpts)::accepted_network_trace
  character*40 :: ibfile
!=======================================================================
!  new variables
!========================================================================
!!  integer, parameter :: nent = 3
  type(ibpt), dimension(:) :: ibary(nent-2)
  type(iapt), dimension(:) :: iXary(nent-2)
  type(cc_augvar), pointer :: cp_lp
  type(cc_bp), pointer :: curr
  integer :: info, llen(nent-2)!, im(2), ict, iloc(4,2), ipl, il
  double precision :: uin(2,nring,nent-2)
!=======================================================================
  PetscErrorCode  ierr
  PetscMPIInt irank, row
  PetscScalar   num
  PetscReal   tol
  PetscBool   setls
  PetscInt   Nin ! # of inner iterations
!=======================================================================
!
#if PETSC_VERSION_LT(3,18,0)
  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
#else
  PetscCallA(PetscInitialize(ierr))
#endif

  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,irank,ierr))

  call InitFSISolve() ! sets up KSP/GMRES + shell matrix for the implicit FSI coupling

  print '("Before reading", 2(i7,1x), 10(e14.6,1x))', ntmax, nfreq, dlt, CFL, nu, bke,bki, vco
  call readpar(CFL,nu,bke,bki,vco,ntmax,nfreq,dlt)
  if(abs(hg-ylength/dble(ny)) > 128.0_dp*epsilon(one)*max(hg,one))then
    error stop 'Cartesian grid requires equal x and y spacing'
  endif
  print '("After reading ", 2(i7,1x), 10(e14.6,1x))', ntmax, nfreq, dlt, CFL, nu, bke,bki, vco
  if(irank.eq.0)then
    write(*,'(a,3(1x,i0),1x,l1,3(1x,es24.16))')'CAMPAIGN_RUNTIME', &
         nx,ny,npts,enforce_fixed_timestep,stage12_initial_network_concentration, &
         stage12_initial_free_concentration,dualchem_pump_width
    write(*,'(a,5(1x,es24.16))')'SIMCELL_PHYSICAL_SCALES', &
         scale_length_um,scale_concentration_millimolar,scale_velocity_um_s, &
         scale_time_s,scale_stress_pa
    write(*,'(a,1x,es24.16)')'SIMCELL_PHYSICAL_INFERRED_VISCOSITY_PA_S', &
         nu*scale_stress_pa*scale_time_s
    if(used_legacy_water_mobility_fallback)then
      write(*,'(a)')'SIMCELL_WARNING legacy kw(1) supplied an omitted water mobility; '// &
           'set water_stress_mobility and water_osmotic_mobility explicitly for a physical map'
    endif
    write(*,'(a,2(1x,es24.16))')'SIMCELL_PHYSICAL_RUN', &
         dlt*scale_time_s,dble(ntmax)*dlt*scale_time_s
    write(*,'(a,8(1x,es24.16))')'SIMCELL_NONDIM_INTERFACE', &
         stage14_adhesion,water_stress_mobility,water_osmotic_mobility, &
         stage12_actin_k_sigma,clstiff,nu,dlt,stage14_external_load_x
  endif

  !!fname='real1.dat'
!  N = nx; M = ny;
  dt = dlt;  ds = two*cpi/dble(nring)
!  dx = hg; dy = hg;

  Nin = 2; Nout=15
!
  print '(1x,"Spacial-x increment: ", 1(e14.6,1x), " Grid nx: ", i5)', hg, nx
  print '(1x,"Spacial-y increment: ", 1(e14.6,1x), " Grid ny: ", i5)', hg, ny
  print '(1x,"viscosity: ", 1(e14.6,1x))', nu
!
  t0 = zero
  tf = .1d0
!  !Nsteps = nint(tf/dt) ! # of time steps determined by dt and tfinal (tf)
!  Nsteps = ntmax! # of time steps determined in input file
!
  print '(1x,"Temporal increment: ", 1(e14.6,1x))', dt
  print '(1x,"Total # of time steps: ", 1(i8,1x))', ntmax
  print *,'============================================================'
  print *,'Starting....'
  print *,'============================================================'
!
!=======================================================================
  call cpu_time(time1) ! initialize time stamp
  time3 = time1
!=======================================================================
! initialize IB locations, and calculate elastic force
  call InitLocCell(xpt, ypt, fsx,fsy, jv) ! (fsx,fsy) stores force, jv store Jacobian factor
  call initialize_dualchem_coupling(xpt,ypt,dualchem_status)
  if (dualchem_status .ne. DUALCHEM_COUPLING_OK) then
    error stop 'Stage 07: dualchem harness initialization failed'
  endif
  call initialize_actin_coefficients(stage12_actin_eta,stage12_actin_eta_s, &
      stage12_actin_k_sigma,stage12_actin_dc,stage12_actin_gamma, &
      stage12_actin_jc,actin_coefficients,actin_status)
  if(actin_status .ne. ACTIN_OK)then
    error stop 'Stage 12: actin coefficient initialization failed'
  endif
  if(irank.eq.0)then
    write(*,'(a,8(1x,es24.16))')'SIMCELL_NONDIM_ACTIN', &
         stage12_actin_eta,stage12_actin_eta_s, &
         actin_coefficients%network_transport_fraction, &
         actin_coefficients%network_diffusivity,stage12_actin_dc, &
         stage12_actin_gamma,stage12_actin_jc,stage12_actin_theta0
    write(*,'(a,5(1x,es24.16))')'SIMCELL_NONDIM_CHEMICAL', &
         dualchem_diffusion,dualchem_kc,dualchem_kp, &
         dualchem_initial_concentration, &
         dualchem_initial_polarization
  endif
  call initialize_actin_oneway(actin_manager,xpt,ypt,actin_coefficients, &
      actin_status)
  if(actin_status .ne. ACTIN_ONEWAY_OK)then
    error stop 'Stage 12: one-way actin manager initialization failed'
  endif
  call get_actin_oneway_mass(actin_manager,actin_mass,actin_status)
  if(actin_status .ne. ACTIN_ONEWAY_OK)then
    error stop 'Stage 12: initial actin mass diagnostic failed'
  endif
  if(irank .eq. 0)then
    write(*,'(a,1x,i8,2(1x,es24.16),1x,i8)') &
        'STAGE12_ACTIN_INITIAL_MASS',0,zero,actin_mass,0
  endif
!
!!!========================================================================
!!! initialize geometry and set initial chemical and actin network 
!!!========================================================================
!!prt  print *, 'step -1'
  !call getGeometry(ibary,llen,-1,uin, dt)
  call getGeometry(ibary,iXary,llen,-1,uin, dt) ! get IB linked list, gridX, normal vector

!!!========================================================================
!
  strlen = len_trim(runname)
!dbg  write(ibfile,'(2a,i4.4)') runname(1:strlen),'.ib.',outcount
!dbg  open(67,file=ibfile,form='formatted',action='write')
!dbg  do i = 1, nring
!dbg    write(67,'(2(e22.14,1x))')xpt(i),ypt(i)
!dbg  enddo
!dbg  write(67,'(2(e22.14,1x))')xpt(1),ypt(1) ! repeat starting point
!dbg  close(67) ! save initial cell configuration
!
  ins = 0 ! for steady stokes 
  select case (ins) ! alpha,beta,gam for setting up stencil. see StokesSolve for details
  case (0) ! RHS NEED TO MULTIPLY dt if dt used in alpha,beta,gam!!!!
    !alpha=4.d0*dt*nu/h/h 
    !beta=dt*nu/h/h
    !gam=dt/h
    alpha=4.d0*nu/hg/hg ! no dt, NEED TO CHECK consistency in RHS 
    beta=nu/hg/hg
    gam=one/hg
  case default
    print *, 'wrong choice for linear solve. stop'
    stop
  end select
  myalpha = alpha; mybeta = beta; mygam = gam
  print '(1x,20(e14.6,1x))', myalpha, mybeta, mygam, alpha, beta, gam
!
! initialize velocity, pressure and source terms
!
  time = zero ! time used in evaluating functions
  cumulative_x_shift=zero
  !!call inituvp(u0, v0, p0, ma1, mb1, uvbc0, 0, dt, time)  ! initialize velocity, uvbc0 is assigned here
  !!call setbc(u0,uvbc0,1)
  !!call setbc(v0,uvbc0,2)
  !!call setbc(p0,uvbc0,3)
  u0 = zero; v0 = zero; p0 = zero; ! set all to initial 0
  uvbc0 = zero
!=======================================================================
! save initial IB locations and u,v, p to data file 
  outcount = 0
  call outuvpx(mpts, u0,v0,p0, xpt,ypt, outcount) ! save initial fluid variables and IB pts
  call get_accepted_physical_chemical_field(physical_chemical_field, &
      dualchem_status)
  if(dualchem_status .ne. DUALCHEM_COUPLING_OK)then
    error stop 'Stage 07: initial physical chemical field unavailable'
  endif
  write(ibfile,'(2a,i4.4)')runname(1:strlen),'.c.',outcount
  open(newunit=snapshot_unit,file=ibfile,access='stream',form='unformatted', &
      status='replace',action='write',iostat=snapshot_iostat)
  if(snapshot_iostat/=0)error stop 'Stage 07: cannot open chemical frame'
  write(snapshot_unit,iostat=snapshot_iostat)physical_chemical_field
  close(snapshot_unit)
  if(snapshot_iostat/=0)error stop 'Stage 07: cannot write chemical frame'
  call get_actin_oneway_fields(actin_manager,physical_network_actin, &
      physical_free_actin,actin_status)
  if(actin_status/=ACTIN_ONEWAY_OK)error stop 'Stage 12: initial fields unavailable'
  if(irank.eq.0)call report_startup(0)
  write(ibfile,'(2a,i4.4)')runname(1:strlen),'.n.',outcount
  open(newunit=snapshot_unit,file=ibfile,access='stream',form='unformatted', &
      status='replace',action='write',iostat=snapshot_iostat)
  if(snapshot_iostat/=0)error stop 'Stage 12: cannot open network frame'
  write(snapshot_unit,iostat=snapshot_iostat)physical_network_actin
  close(snapshot_unit)
  if(snapshot_iostat/=0)error stop 'Stage 12: cannot write network frame'
  write(ibfile,'(2a,i4.4)')runname(1:strlen),'.g.',outcount
  open(newunit=snapshot_unit,file=ibfile,access='stream',form='unformatted', &
      status='replace',action='write',iostat=snapshot_iostat)
  if(snapshot_iostat/=0)error stop 'Stage 12: cannot open free frame'
  write(snapshot_unit,iostat=snapshot_iostat)physical_free_actin
  close(snapshot_unit)
  if(snapshot_iostat/=0)error stop 'Stage 12: cannot write free frame'
!=======================================================================
! Advance in time loop
!=======================================================================
  do it = 1, ntmax
    ! A base interval may be retried as two half intervals.  Only the FSI
    ! trial is speculative; chemistry, actin, geometry, and output are updated
    ! after its marker-normal displacement has passed the one-cell bound.
    accepted_time=dble(it-1)*dt
    remaining_time=dt
    trial_dt=dt
    substep_count=0
    split_count=0
    do while(remaining_time>64.d0*epsilon(one)*max(one,dt))
      time=accepted_time+trial_dt
      accepted_u=u0
      accepted_v=v0
      accepted_p=p0
      accepted_xpt=xpt
      accepted_ypt=ypt
      call getNormal(accepted_xpt,accepted_ypt,nvx,nvy,bdPoly)
    ! Freeze the previously accepted physical trace difference for the entire
    ! FSI solve.  The chemical advance below commits [c]^{n+1} only after this
    ! call returns, so the partitioned order is [c]^n -> FSI^{n+1} -> c^{n+1}.
    call get_accepted_concentration_jump(chemical_jump,dualchem_status)
    if(dualchem_status .ne. DUALCHEM_COUPLING_OK)then
      error stop 'Stage 07: accepted physical concentration jump unavailable'
    end if
    ! Osmotic coupling has its own switch.  It must not be tied to localized
    ! actin polymerization because pump-driven chemical and actin mechanisms
    ! may be enabled simultaneously.
    if(.not.enable_osmotic_feedback)chemical_jump=zero
    ! Stage 14 uses only the previously accepted F-actin state theta_n^n.
    ! The getter returns value copies of the auxiliary field, correction,
    ! physical-side tags, and physical interior trace, so the nonlinear FSI
    ! iteration cannot mutate or observe a partially advanced actin solve.
    call get_actin_oneway_feedback_state(actin_manager, &
         accepted_network_auxiliary,accepted_network_correction, &
         accepted_network_cell_side,accepted_network_trace, &
         actin_source_id,actin_status)
    if(actin_status/=ACTIN_ONEWAY_OK)then
      error stop 'Stage 14: accepted network state unavailable'
    endif
    ! Build the immutable feedback object on the starter MAC layout.  jdu and
    ! jdv are the existing u- and v-face physical-side masks; no new grid,
    ! C++ backend, or interpolation convention is introduced here.
    call build_actin_fsi_feedback(accepted_network_auxiliary, &
         accepted_network_correction,accepted_network_cell_side, &
         jdu(0:nx-1,1:ny),jdv(1:nx,1:ny-1),accepted_network_trace, &
         actin_coefficients,stage12_actin_eta,stage12_actin_k_sigma, &
         actin_source_id,actin_source_id,actin_feedback, &
         actin_feedback_status)
    if(actin_feedback_status/=ACTIN_FEEDBACK_OK)then
      error stop 'Stage 14: active feedback construction failed'
    endif
    ! Morphology acceleration acts only on the spatially varying active stress:
    ! remove its uniform pressure-like part, then amplify the remaining
    ! nonnegative leading-edge stress and its already-gradient-based bulk force.
    if(stage14_shape_force_scale/=one)then
      tmp=minval(actin_feedback%marker_stress)
      actin_feedback%marker_stress=tmp+stage14_shape_force_scale* &
           (actin_feedback%marker_stress-tmp)
      actin_feedback%bulk_force_u=stage14_shape_force_scale* &
           actin_feedback%bulk_force_u
      actin_feedback%bulk_force_v=stage14_shape_force_scale* &
           actin_feedback%bulk_force_v
    endif
    ! Ablate active actin stress and bulk forcing, retaining passive Brinkman
    ! resistance and transported actin. This is an active-driving control.
    if(.not.stage14_active_actin_feedback)then
      actin_feedback%marker_stress=zero
      actin_feedback%bulk_force_u=zero
      actin_feedback%bulk_force_v=zero
    endif
    if(irank.eq.0)write(*,'(a,1x,i0,1x,l1,3(1x,es24.16))') &
      'REDESIGN_ACTIVE_FORCE',it,stage14_active_actin_feedback, &
      maxval(abs(actin_feedback%marker_stress)),maxval(abs(actin_feedback%bulk_force_u)), &
      maxval(abs(actin_feedback%bulk_force_v))
    call AdvanceFSI(u0,v0,p0,xpt,ypt,chemical_jump,actin_feedback,trial_dt)
    ! Use the accepted-time normal, matching the side-change correction's old
    ! geometry.  The absolute value limits motion in either normal direction;
    ! tangential marker redistribution does not force a split.
    max_normal_displacement=maxval(abs((xpt-accepted_xpt)*nvx+ &
         (ypt-accepted_ypt)*nvy))
    max_normal_velocity=max_normal_displacement/trial_dt
    normal_cfl=max_normal_displacement/hg
    if(.not.all(ieee_is_finite(xpt)).or. &
         .not.all(ieee_is_finite(ypt)).or. &
         .not.ieee_is_finite(max_normal_displacement).or. &
         max_normal_displacement>hg)then
      u0=accepted_u
      v0=accepted_v
      p0=accepted_p
      xpt=accepted_xpt
      ypt=accepted_ypt
      if(enforce_fixed_timestep)then
        if(irank.eq.0)write(*,'(a,1x,i0,3(1x,es24.16))') &
             'CAMPAIGN_FIXED_STEP_REJECT',it,time,trial_dt,normal_cfl
        flush(6)
        error stop 'Fixed timestep requires termination instead of subdivision'
      endif
      trial_dt=half*trial_dt
      split_count=split_count+1
      if(split_count>max_normal_step_splits.or. &
           trial_dt<=64.d0*epsilon(one)*max(one,dt))then
        error stop 'Marker-normal timestep subdivision failed to reach CFL <= 1'
      endif
      if(irank.eq.0)write(*,'(a,1x,2(i0,1x),4(es24.16,1x))') &
          'STAGE16_NORMAL_CFL_SPLIT',it,split_count,time, &
          two*trial_dt,max_normal_velocity,normal_cfl
      cycle
    endif
    substep_count=substep_count+1
    if(irank.eq.0)write(*,'(a,1x,2(i0,1x),4(es24.16,1x))') &
        'STAGE16_NORMAL_CFL_ACCEPT',it,substep_count,time,trial_dt, &
        max_normal_velocity,normal_cfl
    call get_fsi_brinkman_diagnostics(brinkman_generation, &
         brinkman_solve_calls,brinkman_inner_iterations,fsi_actin_source_id)
    call get_fsi_solver_diagnostics(fsi_solver_status, &
         fsi_outer_iterations,fsi_initial_residual, &
         fsi_true_absolute_residual,fsi_true_relative_residual, &
         brinkman_krylov_calls,brinkman_worst_reason, &
         brinkman_max_true_absolute_residual, &
         brinkman_max_true_relative_residual)
    ! Optional morphology-only comoving frame.  The x direction is periodic,
    ! so remove accepted rigid translation while retaining every deformation.
    ! Transport sees the corresponding relative velocity; the accumulated lab
    ! displacement remains available as a diagnostic.
    if(stage14_recenter_interface)then
      step_x_shift=sum(xpt)/dble(nring)-xmid
      xpt=xpt-step_x_shift
      cumulative_x_shift=cumulative_x_shift+step_x_shift
      if(irank.eq.0)write(*,'(a,1x,i8,3(1x,es24.16))') &
          'STAGE15_COMOVING_SHIFT',it,time,step_x_shift,cumulative_x_shift
    endif
    if(fsi_actin_source_id/=actin_source_id)then
      error stop 'Stage 14: FSI used a different accepted actin snapshot'
    endif
    max_actin_drag=max(maxval(abs(actin_feedback%drag_u)), &
         maxval(abs(actin_feedback%drag_v)))
    max_actin_bulk_force=max(maxval(abs(actin_feedback%bulk_force_u)), &
         maxval(abs(actin_feedback%bulk_force_v)))
    max_actin_stress=maxval(abs(actin_feedback%marker_stress))
    write(*,'(a,1x,3(i0,1x),3(es24.16,1x),2(i0,1x))') &
         'STAGE14_ACTIN_FEEDBACK',it,actin_source_id,brinkman_generation, &
         max_actin_drag,max_actin_bulk_force,max_actin_stress, &
         brinkman_solve_calls,brinkman_inner_iterations
    write(*,'(a,1x,4(i0,1x))') 'STAGE13_BRINKMAN',it, &
         brinkman_generation,brinkman_solve_calls,brinkman_inner_iterations
    ! AdvanceFSI stores the accepted old markers in the starter-owned
    ! oxpt/oypt arrays before overwriting xpt/ypt.  Chemistry consumes both
    ! time levels and the accepted MAC velocity, but all arguments are INTENT(IN)
    ! in the harness: this call cannot move the interface or feed back to FSI.
    call advance_dualchem_oneway(u0,v0,oxpt,oypt,xpt,ypt,trial_dt,time,dualchem_status)
    if (dualchem_status .ne. DUALCHEM_COUPLING_OK) then
      error stop 'Stage 07: lagged dualchem advance failed'
    endif
    ! The chemical transaction above is the only publisher.  Actin receives a
    ! value copy of that exact accepted FSI velocity/geometry snapshot, then
    ! performs the Stage 11 packed solve.  This publishes A^{m+1} only after
    ! the current FSI solve has finished; the next step may feed A^{m+1} back,
    ! but no current trial/new actin state can enter AdvanceFSI^{m+1}.
    call get_accepted_transport_snapshot(accepted_transport_snapshot, &
        shared_snapshot_id,dualchem_status)
    if(dualchem_status .ne. DUALCHEM_COUPLING_OK)then
      error stop 'Stage 12: accepted shared transport snapshot unavailable'
    endif
    call advance_actin_oneway(actin_manager,accepted_transport_snapshot, &
        shared_snapshot_id,trial_dt,actin_status)
    if(actin_status .ne. ACTIN_ONEWAY_OK)then
      error stop 'Stage 12: one-way actin advance failed'
    endif
    call get_dualchem_mass(dualchem_mass,dualchem_status)
    if (dualchem_status .ne. DUALCHEM_COUPLING_OK) then
      error stop 'Stage 07: dualchem mass diagnostic failed'
    endif
    call get_actin_oneway_mass(actin_manager,actin_mass,actin_status)
    if(actin_status .ne. ACTIN_ONEWAY_OK)then
      error stop 'Stage 12: actin mass diagnostic failed'
    endif
    if(enforce_fixed_timestep)then
      call get_actin_oneway_fields(actin_manager,physical_network_actin, &
           physical_free_actin,actin_status)
      if(actin_status/=ACTIN_ONEWAY_OK)error stop 'Campaign actin fields unavailable'
      call get_accepted_physical_chemical_field(physical_chemical_field,dualchem_status)
      if(dualchem_status/=DUALCHEM_COUPLING_OK)error stop 'Campaign chemical field unavailable'
      if(irank.eq.0)write(*,'(a,1x,i0,5(1x,es24.16))') &
           'CAMPAIGN_COMMITTED',it,time,trial_dt,minval(physical_network_actin), &
           minval(physical_free_actin),minval(physical_chemical_field)
      if(irank.eq.0)call report_startup(it)
      flush(6)
      if(.not.all(ieee_is_finite(physical_network_actin)) .or. &
           .not.all(ieee_is_finite(physical_free_actin)) .or. &
           .not.all(ieee_is_finite(physical_chemical_field))) &
           error stop 'Campaign rejected nonfinite physical concentration'
      if(minval(physical_network_actin)<-1.0e-8_dp .or. &
           minval(physical_free_actin)<-1.0e-8_dp .or. &
           minval(physical_chemical_field)<-1.0e-8_dp) &
           error stop 'Campaign rejected negative physical concentration'
    endif
    ! Re-read only after the complete exterior/interior transaction commits.
    ! These extrema are audit diagnostics for the next step's frozen [c]; the
    ! current FSI solve above used the pre-step copy and cannot see this update.
    call get_accepted_concentration_jump(chemical_jump,dualchem_status)
    if(dualchem_status .ne. DUALCHEM_COUPLING_OK)then
      error stop 'Stage 07: committed physical concentration jump unavailable'
    end if
    if (irank .eq. 0) then
      ! Emit solver diagnostics only after the entire chemical/actin
      ! transaction has committed for this accepted FSI substep.
      write(*,'(a,1x,4(i0,1x),3(es24.16,1x))') 'FIG2_FSI_SOLVER', &
           it,substep_count,fsi_solver_status,fsi_outer_iterations, &
           fsi_initial_residual,fsi_true_absolute_residual, &
           fsi_true_relative_residual
      write(*,'(a,1x,8(i0,1x),2(es24.16,1x))') &
           'FIG2_BRINKMAN_SOLVER',it,substep_count,0, &
           brinkman_generation,brinkman_solve_calls, &
           brinkman_krylov_calls,brinkman_inner_iterations, &
           brinkman_worst_reason,brinkman_max_true_absolute_residual, &
           brinkman_max_true_relative_residual
      write(*,'(a,1x,i8,2(1x,es24.16))') 'STAGE07_DUALCHEM_MASS',it,time,dualchem_mass
      write(*,'(a,1x,i8,3(1x,es24.16))') 'STAGE07_PHYSICAL_JUMP',it,time, &
          minval(chemical_jump),maxval(chemical_jump)
      write(*,'(a,1x,i8,2(1x,es24.16),1x,i8)') &
          'STAGE12_ACTIN_MASS',it,time,actin_mass,shared_snapshot_id
    endif
!========================================================================
! need to update geometry here
!========================================================================
    call getGeometry(ibary,iXary,llen,1,uin,trial_dt)
!
!    call cc_IBbdy(cc_pn, cc_pc, ibary, iXary, llen, uin, 1, time, dt, dif(1))      
!    cc_pc = cc_pn
!
    accepted_time=time
    remaining_time=max(zero,dble(it)*dt-accepted_time)
    if(remaining_time>64.d0*epsilon(one)*max(one,dt))then
      trial_dt=min(trial_dt,remaining_time)
    endif
    enddo
    time=dble(it)*dt

    if (it -nfreq*(it/nfreq) .eq. 0 .and. irank .eq. 0) then ! save data
      !!prt print '(1x, "Step ", i5, " Saving Data ....")', it
      outcount = outcount + 1
      call outuvpx(mpts, u0,v0,p0,xpt,ypt,outcount)
      call get_accepted_physical_chemical_field(physical_chemical_field, &
          dualchem_status)
      if(dualchem_status .ne. DUALCHEM_COUPLING_OK)then
        error stop 'Stage 07: physical chemical output unavailable'
      endif
      write(ibfile,'(2a,i4.4)')runname(1:strlen),'.c.',outcount
      open(newunit=snapshot_unit,file=ibfile,access='stream', &
          form='unformatted',status='replace',action='write', &
          iostat=snapshot_iostat)
      if(snapshot_iostat/=0)error stop 'Stage 07: cannot open chemical frame'
      write(snapshot_unit,iostat=snapshot_iostat)physical_chemical_field
      close(snapshot_unit)
      if(snapshot_iostat/=0)error stop 'Stage 07: cannot write chemical frame'
      call get_actin_oneway_fields(actin_manager,physical_network_actin, &
          physical_free_actin,actin_status)
      if(actin_status/=ACTIN_ONEWAY_OK)error stop 'Stage 12: fields unavailable'
      write(ibfile,'(2a,i4.4)')runname(1:strlen),'.n.',outcount
      open(newunit=snapshot_unit,file=ibfile,access='stream', &
          form='unformatted',status='replace',action='write', &
          iostat=snapshot_iostat)
      if(snapshot_iostat/=0)error stop 'Stage 12: cannot open network frame'
      write(snapshot_unit,iostat=snapshot_iostat)physical_network_actin
      close(snapshot_unit)
      if(snapshot_iostat/=0)error stop 'Stage 12: cannot write network frame'
      write(ibfile,'(2a,i4.4)')runname(1:strlen),'.g.',outcount
      open(newunit=snapshot_unit,file=ibfile,access='stream', &
          form='unformatted',status='replace',action='write', &
          iostat=snapshot_iostat)
      if(snapshot_iostat/=0)error stop 'Stage 12: cannot open free frame'
      write(snapshot_unit,iostat=snapshot_iostat)physical_free_actin
      close(snapshot_unit)
      if(snapshot_iostat/=0)error stop 'Stage 12: cannot write free frame'
    endif
!
    tp3=sum(xpt)/dble(nring); !!tp4 = maxval(xpt);
!    tp1=sum(ypt)/dble(nring)
    tp2 = maxval(xpt)
!    tp4 = xpt(nring/2)
!    call getNormal(xpt,ypt,ndx,ndy,mpoly)
    call curvePoly(bdPoly,tp1,tp4)
    if (irank.eq.0)print '(1x,"It: ", i7, 4(e13.6,1x),20(e14.7,1x))', it, time, tp2, tp3, tp4
!
  enddo ! <-main time evolution loop, it=+1
!=======================================================================
  if(irank .eq. 0)then
    ! Final accepted-state snapshots for the refinement analyzer.  Both files
    ! are header-free binary64 streams in Fortran column-major order:
    ! chemical(i,j), i=1..nx then j=1..ny; jump(k), k=1..mpts.  They are
    ! diagnostics only and are written after all numerical state is accepted.
    call get_accepted_physical_chemical_field(physical_chemical_field, &
        dualchem_status)
    if(dualchem_status .ne. DUALCHEM_COUPLING_OK)then
      error stop 'Stage 07: final physical chemical field unavailable'
    end if
    call get_accepted_concentration_jump(chemical_jump,dualchem_status)
    if(dualchem_status .ne. DUALCHEM_COUPLING_OK)then
      error stop 'Stage 07: final physical concentration jump unavailable'
    end if
    call get_actin_oneway_fields(actin_manager,physical_network_actin, &
        physical_free_actin,actin_status)
    if(actin_status .ne. ACTIN_ONEWAY_OK)then
      error stop 'Stage 12: final physical actin fields unavailable'
    endif
    open(newunit=snapshot_unit,file='./Data/stage07.chemical.final.bin', &
        access='stream',form='unformatted',status='replace',action='write', &
        iostat=snapshot_iostat)
    if(snapshot_iostat/=0)error stop 'Stage 07: cannot open chemical snapshot'
    write(snapshot_unit,iostat=snapshot_iostat)physical_chemical_field
    close(snapshot_unit)
    if(snapshot_iostat/=0)error stop 'Stage 07: cannot write chemical snapshot'
    open(newunit=snapshot_unit,file='./Data/stage07.jump.final.bin', &
        access='stream',form='unformatted',status='replace',action='write', &
        iostat=snapshot_iostat)
    if(snapshot_iostat/=0)error stop 'Stage 07: cannot open jump snapshot'
    write(snapshot_unit,iostat=snapshot_iostat)chemical_jump
    close(snapshot_unit)
    if(snapshot_iostat/=0)error stop 'Stage 07: cannot write jump snapshot'
    ! Header-free binary64, Fortran column-major, matching the Stage 07 field
    ! diagnostics.  Values outside the physical interface are written as zero.
    open(newunit=snapshot_unit,file='./Data/stage12.network.final.bin', &
        access='stream',form='unformatted',status='replace',action='write', &
        iostat=snapshot_iostat)
    if(snapshot_iostat/=0)error stop 'Stage 12: cannot open network snapshot'
    write(snapshot_unit,iostat=snapshot_iostat)physical_network_actin
    close(snapshot_unit)
    if(snapshot_iostat/=0)error stop 'Stage 12: cannot write network snapshot'
    open(newunit=snapshot_unit,file='./Data/stage12.free.final.bin', &
        access='stream',form='unformatted',status='replace',action='write', &
        iostat=snapshot_iostat)
    if(snapshot_iostat/=0)error stop 'Stage 12: cannot open free snapshot'
    write(snapshot_unit,iostat=snapshot_iostat)physical_free_actin
    close(snapshot_unit)
    if(snapshot_iostat/=0)error stop 'Stage 12: cannot write free snapshot'
  end if
!=======================================================================
!
!=======================================================================
! check time spent
!
  if (irank .eq. 0) then
    call cpu_time(time2)
  !!print *, 'Time:', time
    write(*,200)irank, it-1, time2-time3
  endif
!
!dbg  do it = 1, ny
!dbg    print '(1x,1024(i1,1x))', jdu(:,it)
!dbg    !print '(1x,1024(i1,1x))', idu(:,it)
!dbg    !print '(1x,1024(i1,1x))', idf(:,it)
!dbg  enddo
!=======================================================================
  call finalize_dualchem_coupling()
  call FinalizeFSISolve()
  if(irank.eq.0)call report_small_solver_audit()
  if(irank.eq.0)call report_fixed_stencil_audit()
  PetscCallA(PetscFinalize(ierr))
!=======================================================================
200 format(1x,i3, '<<-->>Cputime in ', i7,' step(s):'f16.5,'s')
!=======================================================================
!
contains
  subroutine report_startup(step)
    integer,intent(in)::step
    integer::i,j,k,kp,nf,nr
    real(dp)::cross,area2,cx,cy,cell_area,xx,gf,gr,ff,fr,gt,ft,pg,pf,dx,dy
    area2=0.0_dp;cx=0.0_dp;cy=0.0_dp
    do k=1,nring
      kp=mod(k,nring)+1
      cross=xpt(k)*ypt(kp)-xpt(kp)*ypt(k)
      area2=area2+cross
      cx=cx+(xpt(k)+xpt(kp))*cross
      cy=cy+(ypt(k)+ypt(kp))*cross
    enddo
    cx=cx/(3.0_dp*area2);cy=cy/(3.0_dp*area2)
    cell_area=abs(area2)*0.5_dp*scale_length_um**2
    dx=xlength/real(nx,dp);dy=ylength/real(ny,dp)
    gf=0.0_dp;gr=0.0_dp;ff=0.0_dp;fr=0.0_dp;nf=0;nr=0
    do j=1,ny
      do i=1,nx
        if(abs(physical_network_actin(i,j))+abs(physical_free_actin(i,j))<=0.0_dp)cycle
        xx=xmin+(real(i,dp)-0.5_dp)*dx
        if(xx>=cx)then
          nf=nf+1;gf=gf+physical_free_actin(i,j);ff=ff+physical_network_actin(i,j)
        else
          nr=nr+1;gr=gr+physical_free_actin(i,j);fr=fr+physical_network_actin(i,j)
        endif
      enddo
    enddo
    gt=(gf+gr)*dx*dy;ft=(ff+fr)*dx*dy
    if(min(nf,nr)<=0)error stop 'Startup diagnostic empty actin half'
    gf=gf/real(nf,dp);gr=gr/real(nr,dp);ff=ff/real(nf,dp);fr=fr/real(nr,dp)
    pg=(gf-gr)/max(abs(gf)+abs(gr),epsilon(1.0_dp))
    pf=(ff-fr)/max(abs(ff)+abs(fr),epsilon(1.0_dp))
    ! Actin scale C_A=0.1 mM is fixed throughout this validation campaign.
    write(*,'(a,1x,i0,14(1x,es24.16))') 'F4_STARTUP',step,time*scale_time_s, &
      (cx+cumulative_x_shift)*scale_length_um,cy*scale_length_um,cell_area, &
      ft,gt,ft/max(ft+gt,epsilon(1.0_dp)),0.1_dp*gf,0.1_dp*gr,pg, &
      0.1_dp*ff,0.1_dp*fr,pf,dualchem_pump_start_time*scale_time_s
  end subroutine report_startup

end program IBM
!=======================================================================
