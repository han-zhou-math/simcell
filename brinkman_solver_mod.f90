! Stage 13 variable-coefficient Brinkman inverse.
!
! This module does not introduce a second rectangular Stokes discretization.
! It solves the velocity-space equation
!
!   (I + P0*A_B) u = P0*f,
!
! where P0 is the velocity part of the accepted Fourier/banded incompressible
! Stokes inverse `wrapLinSolve`, and A_B is pointwise nonnegative drag on the
! physical staggered faces.  A final call to wrapLinSolve with f-A_B*u returns
! velocity and pressure in exactly the representation used by the implicit IB
! interpolation and output routines.
!
! The coefficient is copied into module-owned storage before a solve and is
! immutable while the inner PETSc KSP is active.  Stage 14 will construct this
! field from accepted network actin; Stage 13 production freezes zero.
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#ifndef PetscCallA
#define PetscCallA(a) call a
#endif
module brinkman_solver_mod
  use petscsys
  use petscvec
  use petscmat
  use petscksp
  use parameters, only: dp, nx, ny
  use scalar_operator_mod, only: is_finite_scalar
  use linsys, only: wrapLinSolve, myalpha, mybeta, mygam
  implicit none
  private

  integer, parameter, public :: BRINKMAN_OK = 0
  integer, parameter, public :: BRINKMAN_INVALID = 1
  integer, parameter, public :: BRINKMAN_UNINITIALIZED = 2
  integer, parameter, public :: BRINKMAN_SOLVE_FAILED = 3
  integer, parameter, public :: BRINKMAN_NVEL = nx*ny + nx*(ny-1)

  real(dp), save :: frozen_drag_u(nx,ny) = 0.0_dp
  real(dp), save :: frozen_drag_v(nx,ny-1) = 0.0_dp
  integer, save :: frozen_generation = 0
  logical, save :: coefficient_available = .false.
  logical, save :: solver_initialized = .false.
  logical, save :: solve_active = .false.

  type(tKSP), save :: inner_ksp
  type(tMat), save :: inner_mat
  type(tVec), save :: inner_rhs, inner_solution

  public :: initialize_brinkman_solver, finalize_brinkman_solver
  public :: freeze_brinkman_drag, freeze_zero_brinkman_drag
  public :: get_frozen_brinkman_drag, get_brinkman_generation
  public :: pack_physical_velocity, unpack_physical_velocity
  public :: apply_frozen_brinkman_operator, solve_frozen_brinkman

contains

  subroutine initialize_brinkman_solver(status)
    integer, intent(out) :: status
    PetscErrorCode :: ierr
    PetscScalar :: relative_tolerance, absolute_tolerance

    status = BRINKMAN_INVALID
    if (solver_initialized) then
      status = BRINKMAN_OK
      return
    endif

    PetscCallA(KSPCreate(PETSC_COMM_SELF,inner_ksp,ierr))
    PetscCallA(KSPSetType(inner_ksp,KSPGMRES,ierr))
    relative_tolerance = 1.0d-11
    absolute_tolerance = 1.0d-13
    PetscCallA(KSPSetTolerances(inner_ksp,relative_tolerance,absolute_tolerance,PETSC_DEFAULT_REAL,200,ierr))
    PetscCallA(MatCreateShell(PETSC_COMM_SELF,BRINKMAN_NVEL,BRINKMAN_NVEL,BRINKMAN_NVEL,BRINKMAN_NVEL,PETSC_NULL_INTEGER,inner_mat,ierr))
    PetscCallA(MatShellSetOperation(inner_mat,MATOP_MULT,brinkman_mat_mult,ierr))
    PetscCallA(KSPSetOperators(inner_ksp,inner_mat,inner_mat,ierr))
    PetscCallA(VecCreateSeq(PETSC_COMM_SELF,BRINKMAN_NVEL,inner_rhs,ierr))
    PetscCallA(VecDuplicate(inner_rhs,inner_solution,ierr))
    solver_initialized = .true.
    status = BRINKMAN_OK
  end subroutine initialize_brinkman_solver

  subroutine finalize_brinkman_solver(status)
    integer, intent(out) :: status
    PetscErrorCode :: ierr

    status = BRINKMAN_INVALID
    if (solve_active) return
    if (.not. solver_initialized) then
      status = BRINKMAN_OK
      return
    endif
    PetscCallA(VecDestroy(inner_solution,ierr))
    PetscCallA(VecDestroy(inner_rhs,ierr))
    PetscCallA(MatDestroy(inner_mat,ierr))
    PetscCallA(KSPDestroy(inner_ksp,ierr))
    solver_initialized = .false.
    status = BRINKMAN_OK
  end subroutine finalize_brinkman_solver

  subroutine freeze_zero_brinkman_drag(status,generation)
    integer, intent(out) :: status
    integer, intent(out), optional :: generation
    real(dp), allocatable :: drag_u(:,:), drag_v(:,:)

    allocate(drag_u(nx,ny),drag_v(nx,ny-1))
    drag_u = 0.0_dp
    drag_v = 0.0_dp
    call freeze_brinkman_drag(drag_u,drag_v,status,generation)
  end subroutine freeze_zero_brinkman_drag

  subroutine freeze_brinkman_drag(drag_u,drag_v,status,generation)
    real(dp), intent(in) :: drag_u(:,:), drag_v(:,:)
    integer, intent(out) :: status
    integer, intent(out), optional :: generation

    status = BRINKMAN_INVALID
    if (present(generation)) generation = frozen_generation
    if (solve_active) return
    if (size(drag_u,1) /= nx .or. size(drag_u,2) /= ny) return
    if (size(drag_v,1) /= nx .or. size(drag_v,2) /= ny-1) return
    if (.not. all(is_finite_scalar(drag_u))) return
    if (.not. all(is_finite_scalar(drag_v))) return
    if (any(drag_u < 0.0_dp) .or. any(drag_v < 0.0_dp)) return

    frozen_drag_u = drag_u
    frozen_drag_v = drag_v
    coefficient_available = .true.
    frozen_generation = frozen_generation + 1
    if (present(generation)) generation = frozen_generation
    status = BRINKMAN_OK
  end subroutine freeze_brinkman_drag

  subroutine get_frozen_brinkman_drag(drag_u,drag_v,status,generation)
    real(dp), intent(out) :: drag_u(:,:), drag_v(:,:)
    integer, intent(out) :: status
    integer, intent(out), optional :: generation

    status = BRINKMAN_INVALID
    if (present(generation)) generation = frozen_generation
    if (.not. coefficient_available) return
    if (size(drag_u,1) /= nx .or. size(drag_u,2) /= ny) return
    if (size(drag_v,1) /= nx .or. size(drag_v,2) /= ny-1) return
    drag_u = frozen_drag_u
    drag_v = frozen_drag_v
    status = BRINKMAN_OK
  end subroutine get_frozen_brinkman_drag

  subroutine get_brinkman_generation(generation,status)
    integer, intent(out) :: generation,status

    generation = frozen_generation
    if (coefficient_available) then
      status = BRINKMAN_OK
    else
      status = BRINKMAN_UNINITIALIZED
    endif
  end subroutine get_brinkman_generation

  subroutine pack_physical_velocity(u,v,packed,status)
    real(dp), intent(in) :: u(-1:nx+1,-1:ny+1)
    real(dp), intent(in) :: v(-1:nx+1,-1:ny+1)
    real(dp), intent(out) :: packed(:)
    integer, intent(out) :: status
    integer :: i,j,k

    status = BRINKMAN_INVALID
    if (size(packed) /= BRINKMAN_NVEL) return
    if (.not. all(is_finite_scalar(u(0:nx-1,1:ny)))) return
    if (.not. all(is_finite_scalar(v(1:nx,1:ny-1)))) return
    k = 0
    do j=1,ny
      do i=0,nx-1
        k=k+1
        packed(k)=u(i,j)
      enddo
    enddo
    do j=1,ny-1
      do i=1,nx
        k=k+1
        packed(k)=v(i,j)
      enddo
    enddo
    status = BRINKMAN_OK
  end subroutine pack_physical_velocity

  subroutine unpack_physical_velocity(packed,u,v,status)
    real(dp), intent(in) :: packed(:)
    real(dp), intent(inout) :: u(-1:nx+1,-1:ny+1)
    real(dp), intent(inout) :: v(-1:nx+1,-1:ny+1)
    integer, intent(out) :: status
    integer :: i,j,k

    status = BRINKMAN_INVALID
    if (size(packed) /= BRINKMAN_NVEL) return
    if (.not. all(is_finite_scalar(packed))) return
    k = 0
    do j=1,ny
      do i=0,nx-1
        k=k+1
        u(i,j)=packed(k)
      enddo
    enddo
    do j=1,ny-1
      do i=1,nx
        k=k+1
        v(i,j)=packed(k)
      enddo
    enddo
    status = BRINKMAN_OK
  end subroutine unpack_physical_velocity

  subroutine apply_frozen_brinkman_operator(input,output,status)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: output(:)
    integer, intent(out) :: status
    real(dp), allocatable :: u(:,:),v(:,:),force_u(:,:),force_v(:,:)
    real(dp), allocatable :: projected_u(:,:),projected_v(:,:),pressure(:,:)
    real(dp), allocatable :: projected(:)
    integer :: i,j,local_status

    output = 0.0_dp
    status = BRINKMAN_INVALID
    if (.not. coefficient_available) then
      status = BRINKMAN_UNINITIALIZED
      return
    endif
    if (size(input) /= BRINKMAN_NVEL .or. &
        size(output) /= BRINKMAN_NVEL) return
    if (.not. all(is_finite_scalar(input))) return

    allocate(u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1), &
         force_u(-1:nx+1,-1:ny+1),force_v(-1:nx+1,-1:ny+1), &
         projected_u(-1:nx+1,-1:ny+1), &
         projected_v(-1:nx+1,-1:ny+1), &
         pressure(-1:nx+1,-1:ny+1),projected(BRINKMAN_NVEL))

    u = 0.0_dp
    v = 0.0_dp
    call unpack_physical_velocity(input,u,v,local_status)
    if (local_status /= BRINKMAN_OK) return
    force_u = 0.0_dp
    force_v = 0.0_dp
    do j=1,ny
      do i=0,nx-1
        force_u(i,j)=frozen_drag_u(i+1,j)*u(i,j)
      enddo
    enddo
    do j=1,ny-1
      do i=1,nx
        force_v(i,j)=frozen_drag_v(i,j)*v(i,j)
      enddo
    enddo
    call wrapLinSolve(force_u,force_v,projected_u,projected_v,pressure, &
         myalpha,mybeta,mygam)
    call pack_physical_velocity(projected_u,projected_v,projected,local_status)
    if (local_status /= BRINKMAN_OK) return
    output = input + projected
    if (.not. all(is_finite_scalar(output))) then
      output = 0.0_dp
      return
    endif
    status = BRINKMAN_OK
  end subroutine apply_frozen_brinkman_operator

  subroutine solve_frozen_brinkman(f,g,u,v,p,status,iterations,generation_used, &
       krylov_used,convergence_reason,true_absolute_residual, &
       true_relative_residual)
    real(dp), intent(in) :: f(-1:nx+1,-1:ny+1)
    real(dp), intent(in) :: g(-1:nx+1,-1:ny+1)
    real(dp), intent(inout) :: u(-1:nx+1,-1:ny+1)
    real(dp), intent(inout) :: v(-1:nx+1,-1:ny+1)
    real(dp), intent(inout) :: p(-1:nx+1,-1:ny+1)
    integer, intent(out) :: status
    integer, intent(out), optional :: iterations,generation_used
    integer, intent(out), optional :: krylov_used,convergence_reason
    real(dp), intent(out), optional :: true_absolute_residual, &
         true_relative_residual

    real(dp), allocatable :: base_u(:,:),base_v(:,:),base_p(:,:)
    real(dp), allocatable :: trial_u(:,:),trial_v(:,:)
    real(dp), allocatable :: effective_f(:,:),effective_g(:,:)
    real(dp), allocatable :: final_u(:,:),final_v(:,:),final_p(:,:)
    real(dp), allocatable :: packed_rhs(:),packed_solution(:)
    real(dp), allocatable :: operator_image(:)
    PetscErrorCode :: ierr
    PetscInt :: petsc_iterations
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    integer :: reason
#else
    KSPConvergedReason :: reason
#endif
    PetscScalar, pointer :: vector_values(:)
    integer :: i,j,local_status,reason_value
    real(dp) :: absolute_residual,relative_residual,rhs_norm

    status = BRINKMAN_INVALID
    if (present(iterations)) iterations = 0
    if (present(generation_used)) generation_used = frozen_generation
    if (present(krylov_used)) krylov_used = 0
    if (present(convergence_reason)) convergence_reason = 0
    if (present(true_absolute_residual)) true_absolute_residual = 0.0_dp
    if (present(true_relative_residual)) true_relative_residual = 0.0_dp
    if (.not. coefficient_available) then
      status = BRINKMAN_UNINITIALIZED
      return
    endif
    if (.not. all(is_finite_scalar(f(0:nx-1,1:ny)))) return
    if (.not. all(is_finite_scalar(g(1:nx,1:ny-1)))) return

    allocate(base_u(-1:nx+1,-1:ny+1),base_v(-1:nx+1,-1:ny+1), &
         base_p(-1:nx+1,-1:ny+1),trial_u(-1:nx+1,-1:ny+1), &
         trial_v(-1:nx+1,-1:ny+1),effective_f(-1:nx+1,-1:ny+1), &
         effective_g(-1:nx+1,-1:ny+1),final_u(-1:nx+1,-1:ny+1), &
         final_v(-1:nx+1,-1:ny+1),final_p(-1:nx+1,-1:ny+1), &
         packed_rhs(BRINKMAN_NVEL),packed_solution(BRINKMAN_NVEL), &
         operator_image(BRINKMAN_NVEL))

    ! Preserve the accepted direct Stokes path exactly when drag is zero.  This
    ! avoids changing the Stage 12 baseline by a tolerance-sized Krylov error.
    if (all(frozen_drag_u == 0.0_dp) .and. &
        all(frozen_drag_v == 0.0_dp)) then
      call wrapLinSolve(f,g,final_u,final_v,final_p,myalpha,mybeta,mygam)
      if (.not. valid_complete_solution(final_u,final_v,final_p)) return
      u=final_u; v=final_v; p=final_p
      status=BRINKMAN_OK
      return
    endif

    if (.not. solver_initialized) then
      status = BRINKMAN_UNINITIALIZED
      return
    endif

    call wrapLinSolve(f,g,base_u,base_v,base_p,myalpha,mybeta,mygam)
    call pack_physical_velocity(base_u,base_v,packed_rhs,local_status)
    if (local_status /= BRINKMAN_OK) return

#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecGetArrayF90(inner_rhs,vector_values,ierr))
#else
    PetscCallA(VecGetArray(inner_rhs,vector_values,ierr))
#endif
    vector_values(1:BRINKMAN_NVEL)=packed_rhs
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecRestoreArrayF90(inner_rhs,vector_values,ierr))
#else
    PetscCallA(VecRestoreArray(inner_rhs,vector_values,ierr))
#endif
    PetscCallA(VecSet(inner_solution,0.0d0,ierr))
    PetscCallA(KSPSetInitialGuessNonzero(inner_ksp,PETSC_FALSE,ierr))
    solve_active=.true.
    PetscCallA(KSPSolve(inner_ksp,inner_rhs,inner_solution,ierr))
    solve_active=.false.
    PetscCallA(KSPGetConvergedReason(inner_ksp,reason,ierr))
    PetscCallA(KSPGetIterationNumber(inner_ksp,petsc_iterations,ierr))
    if (present(iterations)) iterations=int(petsc_iterations)
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    reason_value=reason
#else
    reason_value=reason%v
#endif
    if(present(krylov_used))krylov_used=1
    if(present(convergence_reason))convergence_reason=reason_value
    if (reason_value <= 0) then
      status=BRINKMAN_SOLVE_FAILED
      return
    endif

#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecGetArrayReadF90(inner_solution,vector_values,ierr))
#else
    PetscCallA(VecGetArrayRead(inner_solution,vector_values,ierr))
#endif
    packed_solution=vector_values(1:BRINKMAN_NVEL)
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecRestoreArrayReadF90(inner_solution,vector_values,ierr))
#else
    PetscCallA(VecRestoreArrayRead(inner_solution,vector_values,ierr))
#endif

    trial_u=0.0_dp
    trial_v=0.0_dp
    call unpack_physical_velocity(packed_solution,trial_u,trial_v,local_status)
    if (local_status /= BRINKMAN_OK) return

    effective_f=f
    effective_g=g
    do j=1,ny
      do i=0,nx-1
        effective_f(i,j)=f(i,j)-frozen_drag_u(i+1,j)*trial_u(i,j)
      enddo
    enddo
    do j=1,ny-1
      do i=1,nx
        effective_g(i,j)=g(i,j)-frozen_drag_v(i,j)*trial_v(i,j)
      enddo
    enddo
    if (.not. all(is_finite_scalar(effective_f(0:nx-1,1:ny)))) return
    if (.not. all(is_finite_scalar(effective_g(1:nx,1:ny-1)))) return
    call wrapLinSolve(effective_f,effective_g,final_u,final_v,final_p, &
         myalpha,mybeta,mygam)
    if (.not. valid_complete_solution(final_u,final_v,final_p)) return
    u=final_u; v=final_v; p=final_p
    status=BRINKMAN_OK

    ! Reapply the frozen matrix only after the returned solution has been
    ! accepted.  This extra operator evaluation is observational: it cannot
    ! change the solution construction or the status returned above.
    absolute_residual=huge(1.0_dp)
    relative_residual=huge(1.0_dp)
    call apply_frozen_brinkman_operator(packed_solution,operator_image, &
         local_status)
    if(local_status==BRINKMAN_OK .and. &
         all(is_finite_scalar(operator_image)))then
      absolute_residual=norm2(operator_image-packed_rhs)
      rhs_norm=norm2(packed_rhs)
      if(rhs_norm>0.0_dp)then
        relative_residual=absolute_residual/rhs_norm
      elseif(absolute_residual==0.0_dp)then
        relative_residual=0.0_dp
      endif
    endif
    if(present(true_absolute_residual)) &
         true_absolute_residual=absolute_residual
    if(present(true_relative_residual)) &
         true_relative_residual=relative_residual
  end subroutine solve_frozen_brinkman

  subroutine brinkman_mat_mult(A,x,y,ierr)
    type(tMat) :: A
    type(tVec) :: x,y
    PetscErrorCode :: ierr
    PetscScalar, pointer :: input_values(:),output_values(:)
    real(dp), allocatable :: input(:),output(:)
    integer :: status

    allocate(input(BRINKMAN_NVEL),output(BRINKMAN_NVEL))
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecGetArrayReadF90(x,input_values,ierr))
#else
    PetscCallA(VecGetArrayRead(x,input_values,ierr))
#endif
    input=input_values(1:BRINKMAN_NVEL)
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecRestoreArrayReadF90(x,input_values,ierr))
#else
    PetscCallA(VecRestoreArrayRead(x,input_values,ierr))
#endif
    call apply_frozen_brinkman_operator(input,output,status)
    if (status /= BRINKMAN_OK) then
      ierr=1
      return
    endif
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecGetArrayF90(y,output_values,ierr))
#else
    PetscCallA(VecGetArray(y,output_values,ierr))
#endif
    output_values(1:BRINKMAN_NVEL)=output
#ifdef SIMCELL_PETSC_LEGACY_ENUM
    PetscCallA(VecRestoreArrayF90(y,output_values,ierr))
#else
    PetscCallA(VecRestoreArray(y,output_values,ierr))
#endif
    ierr=0
  end subroutine brinkman_mat_mult

  logical function valid_complete_solution(u,v,p) result(valid)
    real(dp), intent(in) :: u(-1:nx+1,-1:ny+1)
    real(dp), intent(in) :: v(-1:nx+1,-1:ny+1)
    real(dp), intent(in) :: p(-1:nx+1,-1:ny+1)

    valid=all(is_finite_scalar(u(0:nx,0:ny+1))) .and. &
         all(is_finite_scalar(v(0:nx+1,-1:ny+1))) .and. &
         all(is_finite_scalar(p(0:nx+1,1:ny)))
  end function valid_complete_solution

end module brinkman_solver_mod
