!---------------------------------------------------------------------------
! Implicit fluid-structure interaction (FSI) coupling.
!
! Each timestep solves the implicit-midpoint IB position equation with an
! outer Newton-like iteration.  Its linear correction is solved by PETSc
! GMRES using the matrix-free shell operator zMatMul.  That callback applies
! the implemented approximate Jacobian chain
!
!   elastic-force derivative -> spread -> Stokes solve -> interpolate,
!
! together with the projected water-flux derivative.  Kernel-location,
! normal, metric, and lagged-Brinkman derivatives are frozen, so zMatMul is
! an approximate/semi-Newton Jacobian action rather than the exact derivative
! of every operation in AdvanceFSI.  See the Stage-01 IMPLICIT_IB_METHOD.md
! note for the full notation, derivation, code map, and frozen-term list.
!
! Extracted from fmain.f90's main program (PETSc KSP/MatShell setup and the
! do-ik outer loop) and from zMatMul, which previously lived in linsys.f90.
!---------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscksp.h>
#ifndef PetscCallA
#define PetscCallA(a) call a
#endif
!
module FSIsolve
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use petscksp
  use petscmat
  use petscvec
  use petscsys

  use parameters
  use IBmod
  use IBforce
  use linsys
  use osmotic_feedback_mod, only: OSMOTIC_FEEDBACK_OK, &
      compute_osmotic_slip_velocity
  use brinkman_solver_mod, only: BRINKMAN_OK,initialize_brinkman_solver, &
      finalize_brinkman_solver,solve_frozen_brinkman, &
      get_brinkman_generation,freeze_brinkman_drag
  use fsi_actin_feedback_mod, only: ACTIN_FEEDBACK_OK, &
      actin_fsi_feedback_t,validate_actin_fsi_feedback
  use fsi_actin_force_assembly_mod, only: ACTIN_FORCE_OK, &
      add_actin_fsi_forcing
  use fsi_adhesion_mod, only: ADHESION_OK,compute_adhesion_terms
  use fsi_external_load_mod, only: EXTERNAL_LOAD_OK, &
      compute_external_load_terms

  implicit none

! The input-controlled outer (semi-Newton) iteration count defaults to one: a
! linearly implicit, single-correction update.  Larger values request a more
! fully converged implicit midpoint solve for validation.
! Outer-iteration convergence test on the residual of the midpoint equation
!   F(z) = z - active_step_dt*V(X^n + z/2),   z = X - X^n.
!
! fsrtol is the ACCURACY criterion: the residual must fall fsrtol below its
! ik=1 value. It governs ordinary timesteps, which reach ~1e-11 relative.
!
! fsitol is a SAFE EXIT for when relative reduction is meaningless because the
! problem is already at equilibrium. On a static circle the first residual is
! itself only ~2e-10, and the attainable floor is ~4e-12 -- i.e. 1e-2 RELATIVE,
! which no sane fsrtol would accept, even though 4e-12 absolute is nothing
! (5e-12 per marker, 4e-11 of the radius, ~1e-4 of the temporal discretization
! error). fsitol must therefore sit above the achievable noise floor. Measured
! floors span 1.5e-12..5.3e-11 (validation/10_step1), so 1e-12 was unattainable
! and never fired; 1e-10 clears all of them with margin and is still four
! orders below anything physically resolvable.
  double precision, parameter :: fsitol = 1.d-10 ! absolute: negligible-residual exit
  double precision, parameter :: fsrtol = 1.d-10 ! relative to the ik=1 residual
! Stagnation detection. The attainable residual is limited by a numerical noise
! floor that depends on the marker configuration (it is worst on the exactly
! symmetric analytic initial condition -- see validation/10_step1/RESULTS.md).
! That floor can sit ABOVE fsitol/fsrtol, in which case the tolerance test can
! never be met and the loop burns all Nout iterations in a limit cycle. So also
! stop when the residual has failed to improve on its best value by a factor
! stagfac for Nstag consecutive iterations, and return the best iterate seen.
  double precision, parameter :: stagfac = 0.5d0 ! improvement required to count as progress
  integer, parameter :: Nstag = 3                ! non-improving iterations before giving up
  integer :: nStall = 0                          ! timesteps that ended without converging
  integer :: nStagStep = 0                         ! of those, how many were stagnation

  type(tKSP) :: ksp
  type(tPC) :: pc
  type(tMat) :: mMat
  type(tVec) :: zRHS, xSol
  integer :: active_brinkman_generation = -1
  integer :: step_brinkman_solve_calls = 0
  integer :: step_brinkman_inner_iterations = 0
  integer :: step_brinkman_krylov_calls = 0
  integer :: step_brinkman_worst_reason = 0
  double precision :: step_brinkman_max_true_absolute_residual = 0.d0
  double precision :: step_brinkman_max_true_relative_residual = 0.d0
  integer :: step_fsi_status = 3
  integer :: step_fsi_outer_iterations = 0
  double precision :: step_fsi_initial_residual = 0.d0
  double precision :: step_fsi_true_absolute_residual = 0.d0
  double precision :: step_fsi_true_relative_residual = 0.d0
  type(actin_fsi_feedback_t) :: active_actin_feedback
  integer :: active_actin_source_id = -1
  ! AdvanceFSI sets this before entering the nonlinear solve.  The PETSc shell
  ! callback cannot receive an ordinary Fortran argument, so the residual and
  ! its matrix-free Jacobian share the trial substep size through module state.
  double precision :: active_step_dt = -1.d0

  private
  public :: InitFSISolve, AdvanceFSI, FinalizeFSISolve
  public :: get_fsi_brinkman_diagnostics
  public :: get_fsi_solver_diagnostics

contains
!--------------------------------------------------------------
subroutine InitFSISolve()
  implicit none
  PetscErrorCode :: ierr
  PetscScalar :: num
  integer :: brinkman_status
!
! PETSc's KSP stores the inner Krylov solve used at every outer FSI iterate.
! The linear system has one x and one y unknown per marker (mcoor=2*mpts).
! MatCreateShell allocates no assembled matrix: whenever GMRES needs M*d it
! calls zMatMul, which evaluates the coupled fluid/membrane Jacobian action.
  PetscCallA(KSPCreate(PETSC_COMM_WORLD,ksp,ierr))
  PetscCallA(KSPSetFromOptions(ksp,ierr))
  num = 1.d-8
  PetscCallA(KSPSetTolerances(ksp,num,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,max(200,mcoor),ierr))
  PetscCallA(KSPSetType(ksp,KSPGMRES,ierr))
  ! The shell system has only mcoor marker unknowns.  A full basis avoids the
  ! severe restart stagnation seen for the pressure-scaled JTB coefficients,
  ! while retaining a finite iteration cap for invalid trial geometries.
  PetscCallA(KSPGMRESSetRestart(ksp,mcoor,ierr))
  ! mMat is matrix-free, so factorization preconditioners such as the PETSc
  ! default ILU cannot be constructed.  Keep the outer correction solve
  ! explicitly unpreconditioned on every PETSc build.
  PetscCallA(KSPGetPC(ksp,pc,ierr))
  PetscCallA(PCSetType(pc,PCNONE,ierr))

  PetscCallA(MatCreateShell(PETSC_COMM_SELF,mcoor,mcoor,mcoor,mcoor,PETSC_NULL_INTEGER,mMat,ierr))
  PetscCallA(MatShellSetOperation(mMat,MATOP_MULT,zMatMul,ierr))

  PetscCallA(VecCreateSeq(PETSC_COMM_WORLD,mcoor,zRHS,ierr))
  PetscCallA(VecCreateSeq(PETSC_COMM_WORLD,mcoor,xSol,ierr))
  call initialize_brinkman_solver(brinkman_status)
  if (brinkman_status /= BRINKMAN_OK) then
    error stop 'Stage 13: failed to initialize Brinkman inverse'
  endif
!
  return
end subroutine InitFSISolve
!--------------------------------------------------------------
subroutine AdvanceFSI(u0, v0, p0, xpt, ypt, chemical_jump,actin_feedback,step_dt)
!> Advance the coupled fluid-IB system by one trial timestep of size step_dt,
!> in place: (u0,v0,p0) is overwritten with the new fluid
!> state, (xpt,ypt) with the new boundary position X^{n+1}.
  implicit none
  double precision, dimension(-1:nx+1,-1:ny+1), intent(inout) :: u0, v0, p0
  double precision, dimension(mpts), intent(inout) :: xpt, ypt
  ! Frozen physical trace difference c_i^n-c_e^n.  It is a required input so
  ! the complete implicit IB solve uses exactly one accepted time level.
  double precision, dimension(mpts), intent(in) :: chemical_jump
  type(actin_fsi_feedback_t),intent(in)::actin_feedback
  double precision,intent(in)::step_dt
!
  PetscErrorCode :: ierr
  PetscMPIInt :: row
  PetscScalar :: num
  PetscScalar, pointer :: xx_v(:)
!
  integer :: ik, jk, osmotic_status,actin_feedback_status,actin_force_status
  double precision :: tmp, tp2, tp3, tp4, diagnostic_xamin, diagnostic_yamin
  double precision :: tmp0, tmpbest
  double precision, dimension(mpts) :: xpb, ypb
  integer :: kstag
  logical :: converged, stagnated
  double precision, dimension(-1:nx+1,-1:ny+1) :: u1, v1, ua, va, pa, tiu, tiv
  double precision, dimension(mpts) :: xpa, ypa, vsx,vsy, jv, fsx,fsy, vbx,vby, usx,usy
  double precision, dimension(mpts) :: chemical_slip_x,chemical_slip_y
  double precision, dimension(mpts) :: actin_slip_x,actin_slip_y
  double precision, dimension(mpts) :: adhesion_slip_x,adhesion_slip_y
  integer :: brinkman_status
!
! State-name map used throughout this routine:
!   xpt,ypt = accepted X^n on entry (and X^{n+1} on return)
!   xpa,ypa = current full outer iterate X^{n,k}
!   xpk,ypk = midpoint geometry Xbar^{n,k}=(X^n+X^{n,k})/2
!   u1,v1   = lagged fluid state used by the Brinkman term in velRHS
!   ua,va,pa= Stokes/Brinkman solution on the current midpoint geometry
  if(.not.ieee_is_finite(step_dt).or.step_dt<=zero)then
    error stop 'FSI trial timestep must be positive and finite'
  endif
  active_step_dt=step_dt
  actin_feedback_status=validate_actin_fsi_feedback(actin_feedback)
  if(actin_feedback_status/=ACTIN_FEEDBACK_OK)then
    error stop 'Stage 14: invalid accepted actin feedback object'
  endif
  active_actin_feedback=actin_feedback
  active_actin_source_id=actin_feedback%source_snapshot_id
  call freeze_brinkman_drag(active_actin_feedback%drag_u, &
       active_actin_feedback%drag_v,brinkman_status,active_brinkman_generation)
  if (brinkman_status /= BRINKMAN_OK) then
    error stop 'Stage 14: failed to freeze accepted actin Brinkman drag'
  endif
  step_brinkman_solve_calls=0
  step_brinkman_inner_iterations=0
  step_brinkman_krylov_calls=0
  step_brinkman_worst_reason=0
  step_brinkman_max_true_absolute_residual=0.d0
  step_brinkman_max_true_relative_residual=0.d0
  step_fsi_status=3
  step_fsi_outer_iterations=0
  step_fsi_initial_residual=0.d0
  step_fsi_true_absolute_residual=0.d0
  step_fsi_true_relative_residual=0.d0
  if(.not.all(ieee_is_finite(chemical_jump)))then
    error stop 'Stage 07: nonfinite physical concentration jump'
  end if
  u1 = u0; v1 = v0;
  call getNormal(xpt,ypt,ndx,ndy,bdPoly) ! (ndx,ndy) normal on n level
  call compute_osmotic_slip_velocity(chemical_jump,ndx,ndy, &
      water_osmotic_mobility, &
      chemical_slip_x,chemical_slip_y,osmotic_status)
  if(osmotic_status/=OSMOTIC_FEEDBACK_OK)then
    error stop 'Stage 07: invalid osmotic feedback input'
  end if
  oxpt = xpt; oypt = ypt ! save old IB locs only after input validation
!========================================================================
  xpk = xpt; ypk = ypt ! for \bar X^{n,k} for k= 0
  xpa = xpt; ypa = ypt ! for X^{n,k}
  tmp0 = zero; tmpbest = zero; kstag = 0
  converged = .false.; stagnated = .false.

  do ik = 1, fsi_max_outer_iterations ! Outer iteration for k to reach X^{n+1}=X^{n,k+1}
    step_fsi_outer_iterations=ik
    !(xpk,ypk) store (X^{n,k}+X^n)/2 & X^{n,0} = X^n
    !(xpa,ypa) store X^{n,k}
    ! Evaluate the membrane velocity V(Xbar^{n,k}).  velRHS builds the
    ! Eulerian forcing from the elastic force at Xbar.  solve_fsi_fluid applies
    ! the one frozen coefficient-aware Brinkman inverse (and its exact Stokes
    ! bypass when the coefficient is zero); newInterpS then returns the fluid
    ! velocity at the markers.
    ! RHS fluid part
    call build_fsi_velocity_rhs(u1,v1,xpk,ypk,xpt,ypt,tiu,tiv, &
         actin_slip_x,actin_slip_y,adhesion_slip_x,adhesion_slip_y, &
         actin_force_status)
    if(actin_force_status/=ACTIN_FORCE_OK)then
      error stop 'Stage 14: actin FSI force assembly failed'
    endif
    call solve_fsi_fluid(tiu,tiv,ua,va,pa) ! frozen Brinkman/Stokes inverse
!   One Stokes solve per outer iteration. velRHS reads (u1,v1) only through the
!   Brinkman drag beta(x)*u (bki inside the cell, bke outside), so refreshing it
!   here -- rather than in a second solve after the Newton step -- is the same
!   single drag sweep, one iteration earlier in the lag. With bke=bki=0 it makes
!   no difference at all: velRHS ignores (u1,v1) entirely.
    u1=ua; v1=va
    call newInterpS(xpk,ypk,ua,va,vsx,vsy) ! interplate (ua,va) to IB pts, saved in (vsx,vsy)
    vsx=vsx*hg*hg; vsy=vsy*hg*hg

    ! RHS water flux part
    call IBJacobian(xpk,ypk,jv)
    jv = one/jv !
    call IBForceVec(xpk,ypk,fsx,fsy)
    call normalproject(fsx,fsy, ndx,ndy, vbx,vby)! normal on X^n
    call hadamard(vbx,jv, usx)
    call hadamard(vby,jv, usy) ! RHS on water flux
    ! Because xpk-xpt=(X^{n,k}-X^n)/2, the expression below is
    !
    !   active_step_dt*V(Xbar^{n,k}) - 2*(xpk-xpt)
    !     = active_step_dt*V(Xbar^{n,k}) - (X^{n,k}-X^n)
    !     = -F(z^k),
    !
    ! where z^k=X^{n,k}-X^n and F(z)=z-active_step_dt*V(X^n+z/2).
    ! now set RHS= dt*(fluid+osmosis) - (X^{n,k}-X^n)
    tmp = zero
    do jk = 1, nring
      row = jk-1
      ! With outward n and physical jump [c]=c_i-c_e,
      !
      !   u-X_t = -kw*([c]+Fhat_mem.n)n,
      !   X_t   =  u+kw*Fhat_mem.n*n+kw*[c]*n.
      !
      ! usx/usy already contain the projected elastic force density, while
      ! chemical_slip_* is exactly kw*[c]*n from osmotic_feedback_mod.  The
      ! concentration is a trace value, so unlike the historical pusher force
      ! it must not be multiplied by the inverse interface Jacobian jv.
      num = active_step_dt*(vsx(jk)+water_stress_mobility*usx(jk)+ &
          chemical_slip_x(jk)+actin_slip_x(jk)+adhesion_slip_x(jk))- &
          two*(xpk(jk)-xpt(jk))
      PetscCallA(VecSetValue(zRHS,row,num,INSERT_VALUES,ierr))
      tmp = tmp + (num*num)
      row = jk+nring-1
      num = active_step_dt*(vsy(jk)+water_stress_mobility*usy(jk)+ &
          chemical_slip_y(jk)+actin_slip_y(jk)+adhesion_slip_y(jk))- &
          two*(ypk(jk)-ypt(jk))
      PetscCallA(VecSetValue(zRHS,row,num,INSERT_VALUES,ierr))
      tmp = tmp + (num*num)
    enddo
    tp3=maxval((ua(0:nx-1,1:ny)))
    tp4=maxval((va(1:nx,2:ny-1)))
    tp2=maxval((pa(1:nx,1:ny)))
    tmp = sqrt(tmp)
    print '(1x,i4,1x,"residual",20(e14.6,1x))',ik, tmp, tp3, tp4, tp2
!
!   (xpa,ypa) is the iterate this residual belongs to; remember the best one so a
!   stagnating step can return it instead of an arbitrary point of the cycle.
    if (ik .eq. 1) then
      tmp0 = tmp; tmpbest = tmp; kstag = 0
      xpb = xpa; ypb = ypa
    elseif (tmp < stagfac*tmpbest) then
      tmpbest = tmp; kstag = 0
      xpb = xpa; ypb = ypa
    else
      kstag = kstag + 1
    endif
!   Test convergence HERE, before the Newton step, so we exit holding the very
!   iterate whose residual we just accepted. Testing after the update instead
!   takes one more step, which is harmless while Newton is converging but lands
!   on the bad point of the cycle once the residual is at its noise floor.
    if (tmp < fsitol .or. (ik > 1 .and. tmp < fsrtol*tmp0)) then
      converged = .true.
      exit
    endif
    if (kstag .ge. Nstag) then       ! residual has stopped improving: give up
      stagnated = .true.
      xpa = xpb; ypa = ypb           ! roll back to the best iterate
      do jk = 1, nring
        xpk(jk) = half*(xpa(jk)+xpt(jk))
        ypk(jk) = half*(ypa(jk)+ypt(jk))
      enddo
      call build_fsi_velocity_rhs(u1,v1,xpk,ypk,xpt,ypt,tiu,tiv, &
           actin_slip_x,actin_slip_y,adhesion_slip_x,adhesion_slip_y, &
           actin_force_status)
      if(actin_force_status/=ACTIN_FORCE_OK)then
        error stop 'Stage 14: rollback actin FSI force assembly failed'
      endif
      call solve_fsi_fluid(tiu,tiv,ua,va,pa)
      u1=ua; v1=va
      exit
    endif
    PetscCallA(KSPSetOperators(ksp,mMat,mMat,ierr))
    PetscCallA(KSPSolve(ksp,zRHS,xSol,ierr)) ! sol for system M*(X^{n,k+1}-X^{n,k}) =RHS
!
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecGetArrayReadF90(xSol,xx_v,ierr))
#else
    PetscCallA(VecGetArrayRead(xSol,xx_v,ierr))
#endif
    ! The RHS above is -F(z^k) and mMat is
    ! F'(z) = I - (active_step_dt/2)*J_FSI, so xSol is
    ! the Newton CORRECTION delta = z^{k+1} - z^k, not z^{k+1} itself. It must be
    ! added to the previous iterate (xpa), not to X^n (xpt) -- adding it to xpt
    ! drops the +z^k term, which makes the iteration converge to a point that
    ! does not satisfy F(z)=0. That is invisible whenever the true z is 0 (any
    ! static test) and stalls the residual at ~1e-4 in a moving problem.
    do jk = 1, nring ! update to get X^{n,k+1} = X^{n,k} + delta
      xpa(jk) = xx_v(jk)+xpa(jk); ypa(jk) = xx_v(jk+nring)+ypa(jk)
      xpk(jk) = half*(xpa(jk)+xpt(jk))
      ypk(jk) = half*(ypa(jk)+ypt(jk))
    enddo
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecRestoreArrayReadF90(xSol,xx_v,ierr))
#else
    PetscCallA(VecRestoreArrayRead(xSol,xx_v,ierr))
#endif

  enddo
!
! Both early exits above leave (ua,va,pa) and (u1,v1) consistent with the xpk the
! loop is holding: the convergence test runs before the Newton step, and the
! stagnation branch re-solves after rolling back. Falling out at Nout is the one
! path that does not -- there the loop ended just after a Newton update, so xpk
! has moved past the last solve. Restore consistency for u0/v0/p0.
  if (.not. converged .and. .not. stagnated) then
    call build_fsi_velocity_rhs(u1,v1,xpk,ypk,xpt,ypt,tiu,tiv, &
         actin_slip_x,actin_slip_y,adhesion_slip_x,adhesion_slip_y, &
         actin_force_status)
    if(actin_force_status/=ACTIN_FORCE_OK)then
      error stop 'Stage 14: final actin FSI force assembly failed'
    endif
    call solve_fsi_fluid(tiu,tiv,ua,va,pa)
    u1=ua; v1=va
  endif

  ! Evaluate the nonlinear midpoint residual on the exact state that will be
  ! returned.  In particular this covers the ordinary one-correction path,
  ! whose loop residual belongs to the pre-correction iterate.  These calls
  ! only reconstruct diagnostic vectors and do not participate in acceptance.
  diagnostic_xamin=xamin
  diagnostic_yamin=yamin
  call newInterpS(xpk,ypk,ua,va,vsx,vsy)
  xamin=diagnostic_xamin
  yamin=diagnostic_yamin
  vsx=vsx*hg*hg; vsy=vsy*hg*hg
  call IBJacobian(xpk,ypk,jv)
  jv=one/jv
  call IBForceVec(xpk,ypk,fsx,fsy)
  call normalproject(fsx,fsy,ndx,ndy,vbx,vby)
  call hadamard(vbx,jv,usx)
  call hadamard(vby,jv,usy)
  step_fsi_true_absolute_residual=zero
  do jk=1,nring
    num=active_step_dt*(vsx(jk)+water_stress_mobility*usx(jk)+ &
         chemical_slip_x(jk)+actin_slip_x(jk)+adhesion_slip_x(jk))- &
         two*(xpk(jk)-xpt(jk))
    step_fsi_true_absolute_residual=step_fsi_true_absolute_residual+num*num
    num=active_step_dt*(vsy(jk)+water_stress_mobility*usy(jk)+ &
         chemical_slip_y(jk)+actin_slip_y(jk)+adhesion_slip_y(jk))- &
         two*(ypk(jk)-ypt(jk))
    step_fsi_true_absolute_residual=step_fsi_true_absolute_residual+num*num
  enddo
  step_fsi_true_absolute_residual=sqrt(step_fsi_true_absolute_residual)
  step_fsi_initial_residual=tmp0
  if(tmp0>zero)then
    step_fsi_true_relative_residual=step_fsi_true_absolute_residual/tmp0
  elseif(step_fsi_true_absolute_residual==zero)then
    step_fsi_true_relative_residual=zero
  else
    step_fsi_true_relative_residual=huge(one)
  endif
  if(converged)then
    step_fsi_status=0
  elseif(fsi_max_outer_iterations==1)then
    step_fsi_status=1
  elseif(stagnated)then
    step_fsi_status=2
  else
    step_fsi_status=3
  endif
!
  ! A one-correction run is deliberately a linearly implicit scheme, not a
  ! failed nonlinear solve; retain its state but do not flood its log with
  ! expected non-convergence warnings.
  if (.not. converged .and. fsi_max_outer_iterations>1) then
    nStall = nStall + 1
    if (stagnated) then
      nStagStep = nStagStep + 1
      print '(1x,a,e12.4,a,e12.4,a,i6,a,i6)', &
        & 'WARNING: outer iteration STAGNATED at residual ', tmpbest, &
        & ' (first ', tmp0, '); stagnated/unconverged steps: ', nStagStep, ' /', nStall
    else
      print '(1x,a,i4,a,e12.4,a,e12.4,a,i6)', &
        & 'WARNING: outer iteration did NOT converge in ', fsi_max_outer_iterations, &
        & ' its; residual ', tmp, ' vs first ', tmp0, '; unconverged steps: ', nStall
    endif
  endif
  xpt = (two*xpk-xpt); ypt = (two*ypk-ypt);
  u0 = u1; v0 = v1; p0 = pa
!
  return
end subroutine AdvanceFSI
!--------------------------------------------------------------
subroutine build_fsi_velocity_rhs(u,v,xb,yb,accepted_x,accepted_y,f,g, &
     actin_slip_x,actin_slip_y,adhesion_slip_x,adhesion_slip_y,status)
!> Add the frozen Stage-14 forcing around the unchanged starter velRHS.
!> This helper is called only by the ordinary nonlinear residual.  zMatMul
!> applies the frozen drag through solve_fsi_fluid but deliberately receives
!> none of these constant accepted-state forces.
  implicit none
  double precision,dimension(-1:nx+1,-1:ny+1),intent(in)::u,v
  double precision,dimension(mpts),intent(in)::xb,yb,accepted_x,accepted_y
  double precision,dimension(-1:nx+1,-1:ny+1),intent(out)::f,g
  double precision,dimension(mpts),intent(out)::actin_slip_x,actin_slip_y
  double precision,dimension(mpts),intent(out)::adhesion_slip_x,adhesion_slip_y
  integer,intent(out)::status
  double precision,dimension(mpts)::normal_x,normal_y,jacobian
  double precision,dimension(mpts)::displacement_x,displacement_y
  double precision,dimension(mpts)::adhesion_force_x,adhesion_force_y
  double precision,dimension(mpts)::external_force_x,external_force_y
  double precision,dimension(mpts)::external_slip_x,external_slip_y
  double precision,dimension(-1:nx+1,-1:ny+1)::spread_x,spread_y
  double precision,dimension(7,mpts)::marker_polynomial
  integer::adhesion_status,external_load_status

  status=1
  call velRHS(u,v,xb,yb,f,g)
  call getNormal(xb,yb,normal_x,normal_y,marker_polynomial)
  call IBJacobian(xb,yb,jacobian)
  call add_actin_fsi_forcing(active_actin_feedback,xb,yb,normal_x,normal_y, &
       jacobian,water_stress_mobility,f,g,actin_slip_x,actin_slip_y,status)
  if(status/=ACTIN_FORCE_OK)return
  ! xb is the midpoint (X^n+X^{n,k})/2, hence twice its displacement from
  ! accepted_x is the full implicit marker displacement z^k.
  displacement_x=two*(xb-accepted_x)
  displacement_y=two*(yb-accepted_y)
  call compute_adhesion_terms(displacement_x,displacement_y,normal_x,normal_y, &
       jacobian,active_step_dt,stage14_adhesion,water_stress_mobility, &
       adhesion_force_x,adhesion_force_y,adhesion_slip_x,adhesion_slip_y, &
       adhesion_status)
  if(adhesion_status/=ADHESION_OK)then
    status=1
    return
  endif
  call compute_external_load_terms(normal_x,normal_y,jacobian, &
       two*cpi/dble(nring),stage14_external_load_x,water_stress_mobility, &
       external_force_x,external_force_y,external_slip_x,external_slip_y, &
       external_load_status)
  if(external_load_status/=EXTERNAL_LOAD_OK)then
    status=1
    return
  endif
  adhesion_force_x=adhesion_force_x+external_force_x
  adhesion_force_y=adhesion_force_y+external_force_y
  adhesion_slip_x=adhesion_slip_x+external_slip_x
  adhesion_slip_y=adhesion_slip_y+external_slip_y
  call newSpread(xb,yb,adhesion_force_x,adhesion_force_y,spread_x,spread_y)
  f=f+spread_x
  g=g+spread_y
end subroutine build_fsi_velocity_rhs
!--------------------------------------------------------------
subroutine zMatMul(A,x,y,ierr)
!> Apply the matrix-free approximate Newton matrix M to a marker perturbation.
!>
!> Input x packs d=(dX_x,dX_y).  The callback returns
!>
!>   M d = d - dt I_X L^{-1} S_X[(J_f d)/2 + J_ad d]
!>             - dt[(L_w J_X^{-1}P_nJ_f d)/2 + J_ad,w d],
!>
!> with the code's h^2 interpolation scaling included in I_X.  Here J_f is
!> the elastic-force derivative, S_X/I_X are the IB spread/interpolation
!> operators, L^{-1} is the frozen solve_fsi_fluid inverse, and P_{n^n}
!> projects onto the normal frozen at accepted X^n.  X means the current
!> midpoint xpk/ypk.  For zero drag, L^{-1} reduces exactly to wrapLinSolve.
!>
!> This deliberately omits derivatives of spread/interpolation locations,
!> J_X^{-1}, the normal, and the lagged Brinkman state.  GMRES therefore solves
!> the implemented semi-Newton correction system, not a fully differentiated
!> monolithic residual.  A is unused because PETSc carries the operator through
!> the shell callback itself.
  implicit none
!
  type(tMat) A
  type(tVec) x, y
!
  PetscInt :: row, j
  PetscErrorCode :: ierr
  PetscScalar :: tmp
  PetscScalar, pointer :: xx_v(:)
!
  double precision, dimension(mpts) :: xd, yd, xb,yb, jv, jvx,jvy, vx, vy, xt,yt
  double precision, dimension(mpts) :: normal_x,normal_y,adhesion_force_x, &
       adhesion_force_y,adhesion_slip_x,adhesion_slip_y
  double precision, dimension(7,mpts) :: marker_polynomial
  double precision, dimension(-1:nx+1,-1:ny+1) :: fr, gr, up, vp, pp
  double precision, dimension(-1:nx+1,4) :: uvbc0
  integer :: adhesion_status
!
#ifdef SIMCELL_PETSC_LEGACY_ENUM
  PetscCallA(VecGetArrayReadF90(x,xx_v,ierr)) ! values @ xx_v(1:mcoor)
#else
  PetscCallA(VecGetArrayRead(x,xx_v,ierr)) ! values @ xx_v(1:mcoor)
#endif
!
! PETSc packing convention: all x-marker components, then all y components.
  xd = xx_v(1:mpts)
  yd = xx_v(mpts+1:2*mpts)
  xb = xpk; yb = ypk ! based on X^{n,k}
!
! J_f*d at the midpoint.  IBJacobian supplies the midpoint arclength metric;
! ndx/ndy were computed once at X^n in AdvanceFSI and stay frozen here.
  call IBJacobian(xb,yb, jv) !get Jacobian factor @ IB pts, based on \bar X^{n,k}
  jv = one/jv
  call IBJacobXY(xb,yb, xd,yd, jvx,jvy) ! apply Jacob matrix computed @(xb,yb) to (xd,yd), saved @(jvx,jvy)
!
  ! Semi-Newton choice: project with ndx/ndy from accepted X^n.  The callback
  ! does not differentiate or refresh the normal at the midpoint.
  call normalproject(jvx,jvy, ndx,ndy, xb,yb) !project on normal w/ (xpt,ypt), save to (xb,yb)
  call hadamard(xb,jv, vx) !multiply inverse Jacobian factor to results above
  call hadamard(yb,jv, vy) ! part on water flux, need to multiply kw(1) later

  ! Linearized adhesion uses the full-step perturbation d.  Its reference
  ! force contains J exactly once, while its hydraulic traction contains none.
  call getNormal(xpk,ypk,normal_x,normal_y,marker_polynomial)
  jv=one/jv
  call compute_adhesion_terms(xd,yd,normal_x,normal_y,jv,active_step_dt, &
       stage14_adhesion,water_stress_mobility,adhesion_force_x, &
       adhesion_force_y,adhesion_slip_x,adhesion_slip_y,adhesion_status)
  if(adhesion_status/=ADHESION_OK)error stop 'Implicit adhesion Jacobian failed'
!
! Fluid-mediated derivative: spread J_f*d, solve the homogeneous perturbation
! Stokes/Brinkman problem, then interpolate the result back to the same frozen
! midpoint markers.  A Newton perturbation cannot change prescribed wall
! velocity, so its wall data are exactly homogeneous.  uvbc0 must therefore be
! initialized before the positive-selector setbc calls.
  ! J_f is differentiated with respect to midpoint geometry, giving the 1/2
  ! below.  Adhesion is differentiated with respect to z/dt and has no 1/2.
  jvx=half*jvx+adhesion_force_x
  jvy=half*jvy+adhesion_force_y
  call newSpread(xpk,ypk,jvx,jvy,fr,gr)
  call solve_fsi_fluid(fr,gr,up,vp,pp)
  uvbc0 = zero
  call setbc(up,uvbc0,1)
  call setbc(vp,uvbc0,2)
  call setbc(pp,uvbc0,3)
  call newInterpS(xpk,ypk,up,vp,xt,yt) ! part on fluid velocity
  xt=xt*hg*hg; yt=yt*hg*hg
!
  do j = 1, nring
    tmp = xd(j)-active_step_dt*(xt(j)+half*water_stress_mobility*vx(j)+ &
         adhesion_slip_x(j))
    row = j - 1
    PetscCallA(VecSetValue(y,row,tmp,INSERT_VALUES,ierr))
    tmp = yd(j)-active_step_dt*(yt(j)+half*water_stress_mobility*vy(j)+ &
         adhesion_slip_y(j))
    row = j +nring - 1
    PetscCallA(VecSetValue(y,row,tmp,INSERT_VALUES,ierr))
  enddo
!
#ifdef SIMCELL_PETSC_LEGACY_ENUM
  PetscCallA(VecRestoreArrayReadF90(x,xx_v,ierr))
#else
  PetscCallA(VecRestoreArrayRead(x,xx_v,ierr))
#endif
!
  return
end subroutine zMatMul
!--------------------------------------------------------------
subroutine solve_fsi_fluid(f,g,u,v,p)
!> Apply the one frozen fluid inverse used by both the nonlinear FSI residual
!> and zMatMul.  Keeping this routing in one helper prevents the ordinary and
!> perturbation paths from silently observing different Brinkman fields.
  implicit none
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in) :: f,g
  double precision, dimension(-1:nx+1,-1:ny+1), intent(out) :: u,v,p
  integer :: status,iterations,generation_used,krylov_used,reason
  double precision :: true_absolute_residual,true_relative_residual

  call solve_frozen_brinkman(f,g,u,v,p,status,iterations,generation_used, &
       krylov_used,reason,true_absolute_residual,true_relative_residual)
  if (status /= BRINKMAN_OK) then
    error stop 'Stage 13: frozen Brinkman fluid solve failed'
  endif
  if (generation_used /= active_brinkman_generation) then
    error stop 'Stage 13: Brinkman coefficient changed during FSI solve'
  endif
  step_brinkman_solve_calls=step_brinkman_solve_calls+1
  step_brinkman_inner_iterations=step_brinkman_inner_iterations+iterations
  step_brinkman_krylov_calls=step_brinkman_krylov_calls+krylov_used
  if(krylov_used/=0)then
    if(step_brinkman_worst_reason==0)then
      step_brinkman_worst_reason=reason
    else
      step_brinkman_worst_reason=min(step_brinkman_worst_reason,reason)
    endif
    step_brinkman_max_true_absolute_residual=max( &
         step_brinkman_max_true_absolute_residual,true_absolute_residual)
    step_brinkman_max_true_relative_residual=max( &
         step_brinkman_max_true_relative_residual,true_relative_residual)
  endif
end subroutine solve_fsi_fluid
!--------------------------------------------------------------
subroutine get_fsi_brinkman_diagnostics(generation,solve_calls, &
     inner_iterations,actin_source_id)
  implicit none
  integer, intent(out) :: generation,solve_calls,inner_iterations
  integer, intent(out), optional :: actin_source_id

  generation=active_brinkman_generation
  solve_calls=step_brinkman_solve_calls
  inner_iterations=step_brinkman_inner_iterations
  if(present(actin_source_id))actin_source_id=active_actin_source_id
end subroutine get_fsi_brinkman_diagnostics
!--------------------------------------------------------------
subroutine get_fsi_solver_diagnostics(fsi_status,outer_iterations, &
     initial_residual,true_absolute_residual,true_relative_residual, &
     brinkman_krylov_calls,brinkman_worst_reason, &
     brinkman_max_true_absolute_residual, &
     brinkman_max_true_relative_residual)
  implicit none
  integer,intent(out)::fsi_status,outer_iterations,brinkman_krylov_calls
  integer,intent(out)::brinkman_worst_reason
  double precision,intent(out)::initial_residual,true_absolute_residual
  double precision,intent(out)::true_relative_residual
  double precision,intent(out)::brinkman_max_true_absolute_residual
  double precision,intent(out)::brinkman_max_true_relative_residual

  fsi_status=step_fsi_status
  outer_iterations=step_fsi_outer_iterations
  initial_residual=step_fsi_initial_residual
  true_absolute_residual=step_fsi_true_absolute_residual
  true_relative_residual=step_fsi_true_relative_residual
  brinkman_krylov_calls=step_brinkman_krylov_calls
  brinkman_worst_reason=step_brinkman_worst_reason
  brinkman_max_true_absolute_residual= &
       step_brinkman_max_true_absolute_residual
  brinkman_max_true_relative_residual= &
       step_brinkman_max_true_relative_residual
end subroutine get_fsi_solver_diagnostics
!--------------------------------------------------------------
subroutine FinalizeFSISolve()
  implicit none
  PetscErrorCode :: ierr
  integer :: brinkman_status
!
  call finalize_brinkman_solver(brinkman_status)
  if (brinkman_status /= BRINKMAN_OK) then
    error stop 'Stage 13: failed to finalize Brinkman inverse'
  endif
  PetscCallA(KSPDestroy(ksp,ierr))
!
  return
end subroutine FinalizeFSISolve
!--------------------------------------------------------------
end module FSIsolve
