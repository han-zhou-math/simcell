module linear_solver_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use grid_types
  use advec_diff_solver_mod
  use geometry_mod
  use chemical_interface_flux_mod, only: CHEMICAL_EXTERIOR, CHEMICAL_INTERIOR, &
      CHEMICAL_FLUX_INVALID, CHEMICAL_FLUX_OK, chemical_normalization_t, &
      build_chemical_normalization, build_chemical_robin
  use chemical_pump_profile_mod, only: pnas_gaussian_pump
  use fsi_dualchem_velocity_bridge_mod, only: VELOCITY_BRIDGE_OK, &
      copy_dualchem_face_velocity, get_fsi_velocity_generation
  use parameters, only: dp, nx, ny, npts, xmin, xmax, ymin, ymax, tupi, &
      kc => dualchem_kc, kw => dualchem_kp, scc, dualchem_diffusion, dualchem_pump_start_time, &
      dualchem_interface_gmres_atol_physical, one, zero
  implicit none
  private
  public :: LinearSolver, jmp_func, evaluate_interface_gmres_residual

  integer, parameter, public :: GMRES_NOT_CONVERGED = 0
  integer, parameter, public :: GMRES_RELATIVE_ONLY = 1
  integer, parameter, public :: GMRES_ABSOLUTE_ONLY = 2
  integer, parameter, public :: GMRES_BOTH = 3

  type :: LinearSolver
    ! Solver parameters
    integer :: max_iter
    integer :: restart_m ! Restart parameter (Krylov subspace size)
    real(dp) :: tol
    
    ! Storage for GMRES (Krylov subspace vectors)
    ! v: Basis vectors (n_dof, restart_m+1)
    ! h: Hessenberg matrix (restart_m+1, restart_m)
    real(dp), allocatable :: v(:,:) 
    real(dp), allocatable :: h(:,:) 
    
    ! Auxiliary vectors
    real(dp), allocatable :: c(:), s(:), g(:), y(:)
    
    ! Workspace buffers to avoid repeated allocations
    real(dp), allocatable :: f_grid(:,:), u_grid(:,:)
    real(dp), allocatable :: bc_n(:), bc_s(:), mask_real(:,:)
    real(dp), allocatable :: vx_tmp(:,:), vy_tmp(:,:)
    real(dp), allocatable :: Cxcoef_local(:,:), Cxcoefo(:,:)
    real(dp), allocatable :: x_base(:), y_base(:)  ! Backup for interface positions
    real(dp), allocatable :: phi(:), psi(:), fjmp(:)
    real(dp), allocatable :: vel_ib(:,:)  ! Interface velocity
    real(dp), allocatable :: relative_normal_velocity(:)
    real(dp), allocatable :: eulerian_field(:,:)  ! External field for multigrid solver
    real(dp), allocatable :: u_prev(:,:)  ! Frozen previous step solution (for BDF1 u^n/dt term)
    real(dp) :: dt_step  ! Time step for current solve (frozen during GMRES)
    real(dp), allocatable :: h_zero_work(:)  ! Workspace for caching H(0) during a GMRES solve
    integer(int64) :: kinematics_generation
    logical :: kinematics_ready
    logical :: inner_solve_failed
    logical :: true_residual_ready
    real(dp) :: gmres_initial_normalized_residual
    real(dp) :: gmres_final_normalized_residual
    real(dp) :: gmres_final_relative_residual
    real(dp) :: gmres_diffusion
    real(dp) :: gmres_atol_normalized
    integer :: gmres_convergence_reason
    logical :: gmres_converged
    integer :: gmres_total_krylov_iterations

  contains
    procedure :: init => init_solver
    procedure :: clean => clean_solver
    procedure :: solve => gmres_solve_restart
    procedure :: matvec => matrix_vector_product
    procedure :: compute_H => compute_H_impl
    procedure :: update_history => update_solver_history
    procedure :: prepare_step => prepare_solver_step
    procedure :: set_u_prev => set_previous_step_solution
    procedure :: set_interface_kinematics
    procedure :: get_interface_kinematics
    procedure :: get_last_true_residual
  end type LinearSolver

contains

  pure subroutine evaluate_interface_gmres_residual(final_normalized, &
      initial_normalized,rtol,atol_normalized,converged,reason, &
      relative_residual)
    real(dp),intent(in) :: final_normalized,initial_normalized
    real(dp),intent(in) :: rtol,atol_normalized
    logical,intent(out) :: converged
    integer,intent(out) :: reason
    real(dp),intent(out) :: relative_residual
    logical :: absolute_pass,relative_pass
    real(dp) :: relative_threshold

    converged=.false.
    reason=GMRES_NOT_CONVERGED
    relative_residual=huge(one)
    if(.not.is_finite_backend_scalar(final_normalized) .or. &
       .not.is_finite_backend_scalar(initial_normalized) .or. &
       .not.is_finite_backend_scalar(rtol) .or. &
       .not.is_finite_backend_scalar(atol_normalized))return
    if(final_normalized<zero .or. initial_normalized<zero)return
    if(rtol<=zero .or. atol_normalized<zero)return

    ! rtol<=1 cannot overflow a nonnegative finite initial norm.  For rtol>1,
    ! reject the rounded quotient boundary as well as values above it; the
    ! reviewed boundary test shows that multiplying the quotient itself rounds
    ! to infinity, while its predecessor remains finite.
    if(rtol>one)then
      if(initial_normalized>=huge(one)/rtol)return
    endif
    relative_threshold=rtol*initial_normalized
    if(.not.is_finite_backend_scalar(relative_threshold))return

    if(initial_normalized<=zero)then
      if(final_normalized<=zero)relative_residual=zero
    else
      if(initial_normalized<one)then
        if(final_normalized>huge(one)*initial_normalized)return
      endif
      relative_residual=final_normalized/initial_normalized
      if(.not.is_finite_backend_scalar(relative_residual))return
    endif

    absolute_pass=final_normalized<=atol_normalized
    relative_pass=final_normalized<=relative_threshold
    converged=absolute_pass .or. relative_pass
    if(absolute_pass .and. relative_pass)then
      reason=GMRES_BOTH
    elseif(relative_pass)then
      reason=GMRES_RELATIVE_ONLY
    elseif(absolute_pass)then
      reason=GMRES_ABSOLUTE_ONLY
    endif
  end subroutine evaluate_interface_gmres_residual

  subroutine init_solver(self, max_iter, restart_m, tol, n_dof)
    class(LinearSolver), intent(inout) :: self
    integer, intent(in) :: max_iter, restart_m, n_dof
    real(dp), intent(in) :: tol
    
    self%max_iter = max_iter
    self%restart_m = restart_m
    self%tol = tol
    
    allocate(self%v(n_dof, restart_m + 1))
    allocate(self%h(restart_m + 1, restart_m))
    allocate(self%c(restart_m))
    allocate(self%s(restart_m))
    allocate(self%g(restart_m + 1))
    allocate(self%y(restart_m))
    
    self%v = 0.0_dp
    self%h = 0.0_dp
    self%c = 0.0_dp
    self%s = 0.0_dp
    self%g = 0.0_dp
    self%y = 0.0_dp
    
    ! Allocate workspace buffers
    ! Grid dimensions from parameters module
    allocate(self%f_grid(nx, ny))
    allocate(self%u_grid(nx, ny))
    allocate(self%bc_n(nx))
    allocate(self%bc_s(nx))
    allocate(self%mask_real(nx, ny))
    allocate(self%vx_tmp(nx+1, ny))
    allocate(self%vy_tmp(nx, ny+1))
    allocate(self%Cxcoef_local(npts, 6))
    allocate(self%Cxcoefo(npts, 6))
    allocate(self%x_base(npts))
    allocate(self%y_base(npts))
    allocate(self%phi(npts))
    allocate(self%psi(npts))
    allocate(self%fjmp(npts))
    allocate(self%vel_ib(npts, 2))
    allocate(self%relative_normal_velocity(npts))
    allocate(self%eulerian_field(nx, ny))
    allocate(self%u_prev(nx, ny))
    allocate(self%h_zero_work(n_dof))

    self%f_grid = 0.0_dp
    self%u_grid = 0.0_dp
    self%mask_real = 0.0_dp
    self%vx_tmp = 0.0_dp
    self%vy_tmp = 0.0_dp
    self%Cxcoef_local = 0.0_dp
    self%Cxcoefo = 0.0_dp
    self%x_base = 0.0_dp
    self%y_base = 0.0_dp
    self%phi = 0.0_dp
    self%psi = 0.0_dp
    self%fjmp = 0.0_dp
    self%vel_ib = 0.0_dp
    self%relative_normal_velocity = 0.0_dp
    self%eulerian_field = 0.0_dp
    self%u_prev = 0.0_dp
    self%h_zero_work = 0.0_dp
    self%dt_step = 1.0_dp  ! Default, will be set by caller
    self%kinematics_generation = 0_int64
    self%kinematics_ready = .false.
    self%inner_solve_failed = .false.
    self%true_residual_ready = .false.
    self%gmres_initial_normalized_residual = huge(one)
    self%gmres_final_normalized_residual = huge(one)
    self%gmres_final_relative_residual = huge(one)
    self%gmres_diffusion = huge(one)
    self%gmres_atol_normalized = huge(one)
    self%gmres_convergence_reason = GMRES_NOT_CONVERGED
    self%gmres_converged = .false.
    self%gmres_total_krylov_iterations = 0
  end subroutine init_solver

  subroutine clean_solver(self)
    class(LinearSolver), intent(inout) :: self
    if (allocated(self%v)) deallocate(self%v)
    if (allocated(self%h)) deallocate(self%h)
    if (allocated(self%c)) deallocate(self%c)
    if (allocated(self%s)) deallocate(self%s)
    if (allocated(self%g)) deallocate(self%g)
    if (allocated(self%y)) deallocate(self%y)
    
    ! Deallocate workspace buffers
    if (allocated(self%f_grid)) deallocate(self%f_grid)
    if (allocated(self%u_grid)) deallocate(self%u_grid)
    if (allocated(self%bc_n)) deallocate(self%bc_n)
    if (allocated(self%bc_s)) deallocate(self%bc_s)
    if (allocated(self%mask_real)) deallocate(self%mask_real)
    if (allocated(self%vx_tmp)) deallocate(self%vx_tmp)
    if (allocated(self%vy_tmp)) deallocate(self%vy_tmp)
    if (allocated(self%Cxcoef_local)) deallocate(self%Cxcoef_local)
    if (allocated(self%Cxcoefo)) deallocate(self%Cxcoefo)
    if (allocated(self%x_base)) deallocate(self%x_base)
    if (allocated(self%y_base)) deallocate(self%y_base)
    if (allocated(self%phi)) deallocate(self%phi)
    if (allocated(self%psi)) deallocate(self%psi)
    if (allocated(self%fjmp)) deallocate(self%fjmp)
    if (allocated(self%vel_ib)) deallocate(self%vel_ib)
    if (allocated(self%relative_normal_velocity)) deallocate(self%relative_normal_velocity)
    if (allocated(self%eulerian_field)) deallocate(self%eulerian_field)
    if (allocated(self%u_prev)) deallocate(self%u_prev)
    if (allocated(self%h_zero_work)) deallocate(self%h_zero_work)
    self%kinematics_generation = 0_int64
    self%kinematics_ready = .false.
    self%inner_solve_failed = .false.
    self%true_residual_ready = .false.
    self%gmres_initial_normalized_residual = huge(one)
    self%gmres_final_normalized_residual = huge(one)
    self%gmres_final_relative_residual = huge(one)
    self%gmres_diffusion = huge(one)
    self%gmres_atol_normalized = huge(one)
    self%gmres_convergence_reason = GMRES_NOT_CONVERGED
    self%gmres_converged = .false.
    self%gmres_total_krylov_iterations = 0
  end subroutine clean_solver

  subroutine set_interface_kinematics(self, interface_velocity, r, generation, status)
    class(LinearSolver), intent(inout) :: self
    real(dp), intent(in) :: interface_velocity(:,:), r(:)
    integer(int64), intent(in) :: generation
    integer, intent(out) :: status

    status = CHEMICAL_FLUX_INVALID
    if (.not. allocated(self%vel_ib) .or. &
        .not. allocated(self%relative_normal_velocity)) return
    if (size(interface_velocity,1) /= npts .or. size(interface_velocity,2) /= 2) return
    if (size(r) /= npts) return
    if (.not. all(ieee_is_finite(interface_velocity)) .or. &
        .not. all(ieee_is_finite(r))) return
    if (generation <= 0_int64) return

    self%vel_ib = interface_velocity
    self%relative_normal_velocity = r
    self%kinematics_generation = generation
    self%kinematics_ready = .true.
    status = CHEMICAL_FLUX_OK
  end subroutine set_interface_kinematics

  subroutine get_interface_kinematics(self, interface_velocity, r, generation, status)
    class(LinearSolver), intent(in) :: self
    real(dp), intent(out) :: interface_velocity(:,:), r(:)
    integer(int64), intent(out) :: generation
    integer, intent(out) :: status

    interface_velocity = 0.0_dp
    r = 0.0_dp
    generation = 0_int64
    status = CHEMICAL_FLUX_INVALID
    if (.not. allocated(self%vel_ib) .or. &
        .not. allocated(self%relative_normal_velocity)) return
    if (size(interface_velocity,1) /= npts .or. size(interface_velocity,2) /= 2) return
    if (size(r) /= npts) return
    if (.not. self%kinematics_ready) return

    interface_velocity = self%vel_ib
    r = self%relative_normal_velocity
    generation = self%kinematics_generation
    status = CHEMICAL_FLUX_OK
  end subroutine get_interface_kinematics

  subroutine get_last_true_residual(self, abs_resid, rel_resid, status)
    class(LinearSolver), intent(in) :: self
    real(dp), intent(out) :: abs_resid, rel_resid
    integer, intent(out) :: status

    abs_resid = huge(1.0_dp)
    rel_resid = huge(1.0_dp)
    status = CHEMICAL_FLUX_INVALID
    if (.not. self%true_residual_ready) return
    if (.not. is_finite_backend_scalar( &
        self%gmres_final_normalized_residual) .or. &
        .not. is_finite_backend_scalar( &
        self%gmres_final_relative_residual)) return

    abs_resid = self%gmres_final_normalized_residual
    rel_resid = self%gmres_final_relative_residual
    status = CHEMICAL_FLUX_OK
  end subroutine get_last_true_residual

  ! Restarted GMRES solver
  ! Solves Ax = b where A is the matrix-free operator
  subroutine gmres_solve_restart(self, lag_grid, eul_grid, rhs, solution, kappa, dt, time, isel, first_step, iside_interp, fjmp_in, pde_velocity_zero, steady_state, linear_only, converged_out)
    class(LinearSolver), intent(inout) :: self
    type(LagrangianGrid), intent(inout) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(inout) :: solution(:)
    real(dp), intent(in) :: kappa, dt, time
    integer, intent(in) :: isel
    logical, intent(in) :: first_step
    integer, intent(in) :: iside_interp  ! +1=inside, -1=outside
    real(dp), intent(in), optional :: fjmp_in(:)
    logical, intent(in), optional :: pde_velocity_zero
    logical, intent(in), optional :: steady_state
    logical, intent(in), optional :: linear_only
    ! Stage-06 harness uses this source-compatible trailing result to avoid
    ! committing a partially converged chemical trial.  Legacy callers may
    ! continue to omit it and retain the original printed diagnostics.
    logical, intent(out), optional :: converged_out
    
    integer :: i, j, k, n, iter, m, recursive_reason
    real(dp) :: beta, beta_0, temp, recursive_abs_resid, recursive_rel_resid
    logical :: recursive_cycle_converged,final_predicate_converged
    real(dp), allocatable :: w(:), p_column(:)
    
    n = size(rhs)
    allocate(w(n),p_column(n))
    if (present(converged_out)) converged_out = .false.
    self%inner_solve_failed = .false.
    self%true_residual_ready = .false.
    self%gmres_initial_normalized_residual = huge(one)
    self%gmres_final_normalized_residual = huge(one)
    self%gmres_final_relative_residual = huge(one)
    self%gmres_diffusion = dualchem_diffusion
    self%gmres_atol_normalized = huge(one)
    self%gmres_convergence_reason = GMRES_NOT_CONVERGED
    self%gmres_converged = .false.
    self%gmres_total_krylov_iterations = 0
    recursive_abs_resid = huge(1.0_dp)
    recursive_rel_resid = huge(1.0_dp)
    solution = 0.0_dp ! Initial guess
    m = 0
    beta_0 = -1.0_dp  ! Sentinel: set on first restart

    if(.not.is_finite_backend_scalar(self%gmres_diffusion) .or. &
       self%gmres_diffusion<=zero)goto 900
    if(.not.safe_to_divide_backend( &
       dualchem_interface_gmres_atol_physical,self%gmres_diffusion))goto 900
    self%gmres_atol_normalized=dualchem_interface_gmres_atol_physical/ &
        self%gmres_diffusion
    if(.not.is_finite_backend_scalar(self%gmres_atol_normalized) .or. &
       self%gmres_atol_normalized<zero)goto 900

    ! Compute H(0) once at the start of this GMRES solve.
    ! H(0) is constant throughout the solve because all inputs other than
    ! the Krylov vector p are frozen (geometry, time, dt, boundary conditions, etc.)
    block
      real(dp), allocatable :: zero_vec(:)
      allocate(zero_vec(n))
      zero_vec = 0.0_dp
      call self%compute_H(lag_grid, eul_grid, zero_vec, self%h_zero_work, &
           kappa, dt, time, isel, first_step, iside_interp, fjmp_in, &
           pde_velocity_zero, steady_state)
      deallocate(zero_vec)
    end block

    if (self%inner_solve_failed) goto 900

    ! Outer loop for restarts
    outer_restart: do iter = 1, self%max_iter
        
        ! r0 = b - Ax0
        call self%matvec(lag_grid, eul_grid, solution, w, kappa, dt, time, isel, first_step, iside_interp, fjmp_in, pde_velocity_zero, steady_state, .false.)
        if (self%inner_solve_failed) exit outer_restart
        w = rhs - w
        
        beta = norm2(w)
        ! Capture initial residual norm on first restart for relative tolerance
        if (beta_0 < 0.0_dp) then
            beta_0 = beta
            self%gmres_initial_normalized_residual = beta_0
        endif
        self%gmres_final_normalized_residual = beta
        call evaluate_interface_gmres_residual( &
            self%gmres_final_normalized_residual, &
            self%gmres_initial_normalized_residual,self%tol, &
            self%gmres_atol_normalized,self%gmres_converged, &
            self%gmres_convergence_reason, &
            self%gmres_final_relative_residual)
        self%true_residual_ready = .true.
        
        if (self%gmres_converged) then
            exit outer_restart
        endif
        
        self%v(:, 1) = w / beta
        self%g = 0.0_dp
        self%g(1) = beta
        
        ! Inner Arnoldi loop
        m = self%restart_m
        do j = 1, m
            ! w = A * v_j
            ! GNU Fortran's bounds-checking runtime cannot reliably carry the
            ! descriptor for the allocatable-component section self%v(:,j)
            ! through this polymorphic type-bound call.  Copying the same
            ! Arnoldi column into a contiguous work vector changes neither the
            ! GMRES algebra nor optimized results, and keeps strict Stage-06
            ! validation usable on both macOS and Ascend.
            p_column=self%v(:,j)
            call self%matvec(lag_grid, eul_grid, p_column, w, kappa, dt, time, &
                isel, first_step, iside_interp, fjmp_in, pde_velocity_zero, &
                steady_state, linear_only=.true.)
            if (self%inner_solve_failed) exit outer_restart
            
            ! Arnoldi Process (Modified Gram-Schmidt)
            do i = 1, j
                self%h(i, j) = dot_product(self%v(:, i), w)
                w = w - self%h(i, j) * self%v(:, i)
            end do
            
            self%h(j+1, j) = norm2(w)
            
            if (self%h(j+1, j) > 1.0e-14_dp) then
                self%v(:, j+1) = w / self%h(j+1, j)
            else
                ! Breakdown, but we can still solve
                ! self%h(j+1, j) = 0.0_dp ! Already close to 0
            endif
            
            ! Apply Givens Rotations
            do i = 1, j-1
                temp = self%c(i)*self%h(i, j) + self%s(i)*self%h(i+1, j)
                self%h(i+1, j) = -self%s(i)*self%h(i, j) + self%c(i)*self%h(i+1, j)
                self%h(i, j) = temp
            end do
            
            ! Compute new rotation
            beta = sqrt(self%h(j, j)**2 + self%h(j+1, j)**2)
            if (beta < 1.0e-20_dp) then
                 self%c(j) = 1.0_dp
                 self%s(j) = 0.0_dp
            else
                 self%c(j) = self%h(j, j) / beta
                 self%s(j) = self%h(j+1, j) / beta
            endif
            
            self%h(j, j) = beta
            self%h(j+1, j) = 0.0_dp
            
            ! Update residual vector g
            self%g(j+1) = -self%s(j) * self%g(j)
            self%g(j) = self%c(j) * self%g(j)
            
            recursive_abs_resid = abs(self%g(j+1))
            call evaluate_interface_gmres_residual(recursive_abs_resid, &
                self%gmres_initial_normalized_residual,self%tol, &
                self%gmres_atol_normalized,recursive_cycle_converged, &
                recursive_reason,recursive_rel_resid)
            self%gmres_total_krylov_iterations = &
                self%gmres_total_krylov_iterations+1
            ! print *, "Iter:", (iter-1)*self%restart_m + j, " Rel Resid:", resid_rel
            
            ! The recursive estimate may only shorten this Arnoldi cycle.
            if (recursive_cycle_converged) then
                m = j
                exit
            endif
        end do

        ! Solve upper triangular system Hy = g
        k = m
        self%y(k) = self%g(k) / self%h(k, k)
        do i = k-1, 1, -1
            self%y(i) = (self%g(i) - dot_product(self%h(i, i+1:k), self%y(i+1:k))) / self%h(i, i)
        end do
        
        ! Update solution: x = x + V_k * y
        do i = 1, k
            solution = solution + self%y(i) * self%v(:, i)
        end do

        ! The recursive GMRES estimate is not an acceptance criterion.  Apply
        ! the matrix-free operator to the updated candidate and store the true
        ! residual used by the caller's transaction gate.
        call self%matvec(lag_grid, eul_grid, solution, w, kappa, dt, time, &
            isel, first_step, iside_interp, fjmp_in, pde_velocity_zero, &
            steady_state, linear_only=.true.)
        if (self%inner_solve_failed) exit outer_restart
        self%gmres_final_normalized_residual = norm2(rhs-w)
        call evaluate_interface_gmres_residual( &
            self%gmres_final_normalized_residual, &
            self%gmres_initial_normalized_residual,self%tol, &
            self%gmres_atol_normalized,self%gmres_converged, &
            self%gmres_convergence_reason, &
            self%gmres_final_relative_residual)
        self%true_residual_ready = .true.
        if (self%gmres_converged) exit outer_restart
        
    end do outer_restart
    
900 continue
    ! Re-evaluate the stored true residual at the final return boundary.  The
    ! recursive estimate is deliberately absent from this acceptance path.
    call evaluate_interface_gmres_residual( &
        self%gmres_final_normalized_residual, &
        self%gmres_initial_normalized_residual,self%tol, &
        self%gmres_atol_normalized,final_predicate_converged, &
        self%gmres_convergence_reason,self%gmres_final_relative_residual)
    self%gmres_converged=final_predicate_converged .and. &
        self%true_residual_ready .and. .not.self%inner_solve_failed
    if(.not.self%gmres_converged) &
        self%gmres_convergence_reason=GMRES_NOT_CONVERGED
    if (.not.self%gmres_converged) then
        print *, '***********************************************************'
        print *, 'WARNING: GMRES did not converge to tolerance'
        print *, '  Target tolerance (relative): ', self%tol
        print *, '  Target tolerance (absolute): ', self%gmres_atol_normalized
        print *, '  Initial residual:            ', self%gmres_initial_normalized_residual
        print *, '  Final true absolute residual:', self%gmres_final_normalized_residual
        print *, '  Final true relative residual:', self%gmres_final_relative_residual
        print *, '  Inner solve failed:          ', self%inner_solve_failed
        print *, '  Iterations used:             ', self%gmres_total_krylov_iterations
        print *, '***********************************************************'
    else
        print *, 'GMRES converged. True rel resid: ', &
            self%gmres_final_relative_residual, &
            ' True abs resid: ', self%gmres_final_normalized_residual, &
            ' reason: ', self%gmres_convergence_reason, &
            ' iterations: ', self%gmres_total_krylov_iterations
        if (present(converged_out)) converged_out = .true.
    endif
    
    deallocate(w,p_column)

  end subroutine gmres_solve_restart

  ! Compute H(psi): The full affine operator
  ! H(psi) = A_linear(psi) + constant_term
  subroutine compute_H_impl(self, lag_grid, eul_grid, psi, h_out, kappa, dt, time, isel, first_step, iside_interp, fjmp_in, pde_velocity_zero, steady_state, u_at_interface_in, u_at_interface_out)
    class(LinearSolver), intent(inout) :: self
    type(LagrangianGrid), intent(inout) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    real(dp), intent(in) :: psi(:) 
    real(dp), intent(out) :: h_out(:)
    real(dp), intent(in) :: kappa, dt, time
    integer, intent(in) :: isel
    logical, intent(in) :: first_step
    integer, intent(in) :: iside_interp
    real(dp), intent(in) :: fjmp_in(:)
    logical, intent(in), optional :: pde_velocity_zero
    logical, intent(in), optional :: steady_state
    ! Optional output: one-sided field values at interface (for coupling optimization)
    real(dp), intent(out), optional :: u_at_interface_in(:)   ! Interior-side values (iside=+1)
    real(dp), intent(out), optional :: u_at_interface_out(:)  ! Exterior-side values (iside=-1)

    
    integer :: iter_cpp, i, bridge_status, flux_status,correction_status
    integer(int64) :: bridge_generation
    integer :: nx_local, ny_local
    real(dp) :: avec(6), tp1, tp2, kappa_solve, kappa_normalized
    real(dp) :: inverse_dt_effective
    real(dp) :: jmp_coef, coef_u, coef_grad, opposite_coef, shift_constant
    real(dp) :: val_in, val_out  ! Store one-sided values temporarily
    real(dp) :: source_jump_normalized(size(fjmp_in))
    type(chemical_normalization_t) :: normalization
    logical :: inner_converged,zero_velocity_requested,steady_state_requested

    nx_local = eul_grid%nx_grid
    ny_local = eul_grid%ny_grid
    h_out = 0.0_dp
    if (present(u_at_interface_in)) u_at_interface_in = 0.0_dp
    if (present(u_at_interface_out)) u_at_interface_out = 0.0_dp
    zero_velocity_requested=.false.
    if(present(pde_velocity_zero))zero_velocity_requested=pde_velocity_zero
    steady_state_requested=.true.
    if(present(steady_state))steady_state_requested=steady_state

    ! The backend is unit-diffusion.  Normalize every disposable Cartesian
    ! input from the one accepted physical D; physical dt, time, geometry, and
    ! the frozen FSI snapshot remain unchanged at the coupling boundary.
    call build_chemical_normalization(dualchem_diffusion,dt,normalization, &
        flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK .or. &
        .not.is_finite_backend_scalar(time))then
        self%inner_solve_failed=.true.
        return
    endif
    if(.not.safe_to_scale_backend(kappa,normalization%inverse_diffusion) .or. &
        .not.all(safe_to_scale_backend(fjmp_in,normalization%inverse_diffusion)) .or. &
        .not.safe_to_divide_backend(one,normalization%dt_effective))then
        self%inner_solve_failed=.true.
        return
    endif
    kappa_normalized=kappa*normalization%inverse_diffusion
    source_jump_normalized=fjmp_in*normalization%inverse_diffusion
    inverse_dt_effective=one/normalization%dt_effective

    ! Every affine-operator evaluation must consume the exact generation
    ! frozen by the coupling transaction.  A stale bridge is an inner failure,
    ! never a process-aborting error or a usable Krylov vector.
    if (.not. self%kinematics_ready) then
        self%inner_solve_failed = .true.
        return
    endif

    call get_fsi_velocity_generation(bridge_generation, bridge_status)
    if (bridge_status /= VELOCITY_BRIDGE_OK .or. &
        bridge_generation /= self%kinematics_generation) then
        self%inner_solve_failed = .true.
        return
    endif

    ! 0. CRITICAL: Recompute correction coefficients for this psi value
    !    This is the key step that makes the matrix-free matvec work correctly.
    !    Without this, the PDE solution doesn't depend on psi, so A*p ≈ 0.5*p.
    !    Note: fjmp_in must be used instead of self%fjmp since prepare_step may not be called.
    
    ! DEBUG MODE (commented out - testing Option 3):
    ! Analytical correction: crc(i,j) = u_in(x,y) - u_out(x,y) = [u] at irregular points
!    block
!        integer :: ii, jj
!        real(dp) :: xx, yy, u_in_val, u_out_val
!        eul_grid%crc = 0.0_dp
!        do jj = 1, ny_local
!            yy = eul_grid%y_min + (real(jj, dp) - 0.5_dp) * eul_grid%dy
!            do ii = 1, nx_local
!                xx = eul_grid%x_min + (real(ii, dp) - 0.5_dp) * eul_grid%dx
!                if (eul_grid%id(ii, jj) == 1 .or. eul_grid%id(ii, jj) == -1) then
!                    u_in_val = cos(tupi * xx) * sin(tupi * yy)  ! get_u_exact
!                    u_out_val = cos(tupi * xx) * sin(yy)         ! get_u_out
!                    eul_grid%crc(ii, jj) = u_in_val - u_out_val
!                endif
!            enddo
!        enddo
!    end block
    
    ! NORMAL MODE: Use getCorrection to compute correction coefficients
!    if (present(fjmp_in)) then
        if(zero_velocity_requested)then
            call getCorrection(lag_grid,lag_grid%phi,psi, &
                source_jump_normalized,normalization%dt_effective,time, &
                self%Cxcoef_local,self%Cxcoefo,kappa_normalized, &
                pde_velocity_zero=.true.,steady_state=steady_state_requested, &
                status=correction_status)
        else
            call getCorrection(lag_grid,lag_grid%phi,psi, &
                source_jump_normalized,normalization%dt_effective,time, &
                self%Cxcoef_local,self%Cxcoefo,kappa_normalized, &
                steady_state=steady_state_requested, &
                pde_velocity_scale=normalization%inverse_diffusion, &
                status=correction_status)
        endif
        if(correction_status/=CHEMICAL_FLUX_OK)then
            self%inner_solve_failed=.true.
            return
        endif
!    else
!        call getCorrection(lag_grid, lag_grid%phi, psi, self%fjmp, dt, time, &
!                          self%Cxcoef_local, self%Cxcoefo, kappa, pde_velocity_zero, steady_state)
!    endif
    call apply_full_correction(lag_grid, eul_grid, self%Cxcoef_local, &
        self%Cxcoefo,iside_interp,steady_state_requested, &
        normalization%dt_effective)
    
    ! 1. Set source from external Eulerian field + u_prev/dt for time-dependent
    !    u_prev is FROZEN during all GMRES iterations within a time step
    if(.not.all(safe_to_scale_backend(self%eulerian_field, &
        normalization%inverse_diffusion)))then
        self%inner_solve_failed=.true.
        return
    endif
    self%f_grid=self%eulerian_field*normalization%inverse_diffusion
    if (allocated(self%u_prev) .and. .not.steady_state_requested) then
        if(.not.all(safe_to_scale_backend(self%u_prev,inverse_dt_effective)))then
            self%inner_solve_failed=.true.
            return
        endif
        self%u_grid=self%u_prev*inverse_dt_effective
        if(.not.all(safe_to_add_backend(self%f_grid,self%u_grid)))then
            self%inner_solve_failed=.true.
            return
        endif
        self%f_grid=self%f_grid+self%u_grid
    endif
    
    ! 2. Apply freshly cleared correction to f_grid (RHS)
    !    This accounts for points that changed status due to interface movement
    if (allocated(self%Cxcoefo)) then
        call apply_fresh_cleared_correction(lag_grid, eul_grid, self%Cxcoefo, self%f_grid, &
            iside_interp,steady_state_requested,normalization%dt_effective)
        if(.not.all(ieee_is_finite(self%f_grid)))then
            self%inner_solve_failed=.true.
            return
        endif
    endif
    
    eul_grid%u = 0.0_dp ! initial guess for u
    
    ! 2. Solve Eulerian PDE: -Δu + v\cdot \nabla u = f_eulerian + corrections
    self%u_grid = 0.0_dp
    ! NOTE: bc_n and bc_s should be set by the caller (test program) - do not reset here
    ! self%bc_n = 0.0_dp
    ! self%bc_s = 0.0_dp
    ! C++ interface expects a 0/1 mask (1 = interior), so convert idf (+/-1) accordingly.
    self%mask_real = 0.5_dp * (real(eul_grid%idf(1:nx_local, 1:ny_local), dp) + 1.0_dp)
    
    if (zero_velocity_requested) then
        self%vx_tmp = 0.0_dp
        self%vy_tmp = 0.0_dp
    else ! use the same frozen FSI state at all Cartesian faces
        call copy_dualchem_face_velocity(self%vx_tmp, self%vy_tmp, bridge_status)
        if (bridge_status /= VELOCITY_BRIDGE_OK) then
            self%inner_solve_failed = .true.
            return
        endif
        if(.not.all(safe_to_scale_backend(self%vx_tmp, &
            normalization%inverse_diffusion)) .or. &
            .not.all(safe_to_scale_backend(self%vy_tmp, &
            normalization%inverse_diffusion)))then
            self%inner_solve_failed=.true.
            return
        endif
        self%vx_tmp=self%vx_tmp*normalization%inverse_diffusion
        self%vy_tmp=self%vy_tmp*normalization%inverse_diffusion
    endif
    
    ! 3. Solve advection-diffusion equation with correction terms
    if (.not.steady_state_requested) then
        if(.not.safe_to_add_backend(kappa_normalized,inverse_dt_effective))then
            self%inner_solve_failed=.true.
            return
        endif
        kappa_solve=kappa_normalized+inverse_dt_effective
    else
        kappa_solve=kappa_normalized
    endif
    if(.not.ieee_is_finite(kappa_solve))then
        self%inner_solve_failed=.true.
        return
    endif
    call solve_advec_diff(xmin, ymin, xmax, ymax, nx_local, ny_local, &
                          self%vx_tmp, self%vy_tmp, &
                          kappa_solve, self%f_grid, self%bc_n, self%bc_s, &
                          self%mask_real, eul_grid%crc, &
                          self%u_grid, iter_cpp, inner_converged)
    if (.not. inner_converged) then
        self%inner_solve_failed = .true.
        return
    endif
    
    ! 5. Interpolate Eulerian solution to Lagrangian grid for averaged jump in value and normal gradient
    eul_grid%u(1:nx_local, 1:ny_local) = self%u_grid 

    do i = 1, npts
        call valatIBpt(avec, i, self%Cxcoef_local, lag_grid, eul_grid, eul_grid%u, 1)
        val_in = avec(1)
        tp1 = avec(1)
        tp2 = lag_grid%normal(i,1)*avec(2) + lag_grid%normal(i,2)*avec(3)
        
        call valatIBpt(avec, i, self%Cxcoef_local, lag_grid, eul_grid, eul_grid%u, -1)
        val_out = avec(1)
        tp1 = tp1 + avec(1)
        tp2 = tp2 + lag_grid%normal(i,1)*avec(2) + lag_grid%normal(i,2)*avec(3)
        
        ! Save one-sided values if requested
        if (present(u_at_interface_in)) u_at_interface_in(i) = val_in
        if (present(u_at_interface_out)) u_at_interface_out(i) = val_out
        
        ! Compute jmp coefficient (active pump term)
        jmp_coef = jmp_func(scc(i), time)
        
        call build_chemical_robin(iside_interp, &
            self%relative_normal_velocity(i),kc,jmp_coef, &
            normalization%diffusion, &
            coef_u,coef_grad,opposite_coef,shift_constant,flux_status)
        if (flux_status /= CHEMICAL_FLUX_OK) then
            self%inner_solve_failed = .true.
            h_out = 0.0_dp
            return
        endif

        ! Apply the corrected moving-interface Robin coefficients (explicitly
        ! diffusion-normalized) to the existing average-plus-jump reconstruction.
        ! The opposite trace and stored-variable shift are assembled once by the
        ! harness in lag_grid%rhs.
        if (iside_interp == CHEMICAL_INTERIOR) then
            h_out(i) = coef_u * (0.5_dp * tp1) + coef_grad * (0.5_dp * tp2) - 0.5_dp*psi(i) - lag_grid%rhs(i) 
            ! DEBUG: Print for first interface point
            !if (i == 1) then
            !    print *, "DEBUG INTERIOR i=1: tp1=", tp1, "tp2=", tp2
            !    print *, "  coef_u=", coef_u, "psi=", psi(i), "rhs=", lag_grid%rhs(i)
            !    print *, "  h_out=", h_out(i)
            !endif
        else if (iside_interp == CHEMICAL_EXTERIOR) then
            h_out(i) = coef_u * (0.5_dp * tp1) + coef_grad * (0.5_dp * tp2) + 0.5_dp*psi(i) - lag_grid%rhs(i)
            ! DEBUG: Print for first interface point
            !if (i == 1) then
            !    print *, "DEBUG EXTERIOR i=1: tp1=", tp1, "tp2=", tp2
            !    print *, "  coef_u=", coef_u, "psi=", psi(i), "rhs=", lag_grid%rhs(i)
            !    print *, "  h_out=", h_out(i)
            !endif
        else
            self%inner_solve_failed = .true.
            h_out = 0.0_dp
            return
        endif
    enddo
    
  end subroutine compute_H_impl

  ! Matrix-Free Matrix-Vector Product: w = A * p
  ! Constructed using H(psi): w = H(p) - H(0)
  ! H(0) is precomputed once at the start of gmres_solve_restart and stored
  ! in self%h_zero_work, since all inputs other than p are frozen during the solve.
  subroutine matrix_vector_product(self, lag_grid, eul_grid, p, w, kappa, dt, time, isel, first_step, iside_interp, fjmp_in, pde_velocity_zero, steady_state, linear_only)
    class(LinearSolver), intent(inout) :: self
    type(LagrangianGrid), intent(inout) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    real(dp), intent(in) :: p(:) ! Input vector (on Lagrangian grid)
    real(dp), intent(out) :: w(:) ! Output vector (on Lagrangian grid)
    real(dp), intent(in) :: kappa, dt, time
    integer, intent(in) :: isel
    logical, intent(in) :: first_step
    integer, intent(in) :: iside_interp  ! +1=inside, -1=outside
    real(dp), intent(in), optional :: fjmp_in(:)
    logical, intent(in), optional :: pde_velocity_zero
    logical, intent(in), optional :: steady_state
    logical, intent(in), optional :: linear_only

    real(dp), allocatable :: h_psi(:)

    allocate(h_psi(size(p)))

    ! Compute H(p)
    call self%compute_H(lag_grid, eul_grid, p, h_psi, kappa, dt, time, isel, first_step, iside_interp, fjmp_in, pde_velocity_zero, steady_state)

    ! w = H(p) - H(0), using precomputed H(0) from self%h_zero_work
    w = h_psi - self%h_zero_work

    deallocate(h_psi)

  end subroutine matrix_vector_product

  subroutine update_solver_history(self)
    class(LinearSolver), intent(inout) :: self
    
    if (allocated(self%Cxcoef_local) .and. allocated(self%Cxcoefo)) then
        self%Cxcoefo = self%Cxcoef_local
    endif
    
  end subroutine update_solver_history

  subroutine prepare_solver_step(self, lag_grid, eul_grid, psi, kappa, dt, time, fjmp_in, pde_velocity_zero, steady_state, phi_in, iside)
    class(LinearSolver), intent(inout) :: self
    type(LagrangianGrid), intent(inout) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    real(dp), intent(in) :: psi(:) 
    real(dp), intent(in) :: kappa, dt, time
    real(dp), intent(in), optional :: fjmp_in(:)
    logical, intent(in), optional :: pde_velocity_zero
    logical, intent(in), optional :: steady_state
    real(dp), intent(in), optional :: phi_in(:)
    integer, intent(in), optional :: iside
    type(chemical_normalization_t)::normalization
    real(dp)::kappa_normalized
    integer::flux_status,correction_status
    logical::zero_velocity_requested,steady_state_requested
    
    ! 1. Setup Data
     if (present(phi_in)) then
        self%phi = phi_in
    elseif (allocated(lag_grid%phi)) then
        self%phi = lag_grid%phi
    else
        self%phi = 0.0_dp
    endif
    
    self%psi = psi
    
    if (present(fjmp_in)) then
        self%fjmp = fjmp_in
    else
        self%fjmp = 0.0_dp
    endif

    zero_velocity_requested=.false.
    if(present(pde_velocity_zero))zero_velocity_requested=pde_velocity_zero
    steady_state_requested=.true.
    if(present(steady_state))steady_state_requested=steady_state
    call build_chemical_normalization(dualchem_diffusion,dt,normalization,flux_status)
    if(flux_status/=CHEMICAL_FLUX_OK .or. &
        .not.is_finite_backend_scalar(time) .or. &
        .not.safe_to_scale_backend(kappa,normalization%inverse_diffusion) .or. &
        .not.all(safe_to_scale_backend(self%fjmp,normalization%inverse_diffusion)))then
      self%inner_solve_failed=.true.
      return
    endif
    kappa_normalized=kappa*normalization%inverse_diffusion
    self%fjmp=self%fjmp*normalization%inverse_diffusion
    
    ! 2. Compute Correction Coefficients (Geometric + Solution Dependent)
    if(zero_velocity_requested)then
      call getCorrection(lag_grid,self%phi,self%psi,self%fjmp, &
          normalization%dt_effective,time,self%Cxcoef_local,self%Cxcoefo, &
          kappa_normalized,pde_velocity_zero=.true., &
          steady_state=steady_state_requested,status=correction_status)
    else
      call getCorrection(lag_grid,self%phi,self%psi,self%fjmp, &
          normalization%dt_effective,time,self%Cxcoef_local,self%Cxcoefo, &
          kappa_normalized,steady_state=steady_state_requested, &
          pde_velocity_scale=normalization%inverse_diffusion, &
          status=correction_status)
    endif
    if(correction_status/=CHEMICAL_FLUX_OK)then
      self%inner_solve_failed=.true.
      return
    endif
                      
    ! 3. Apply Correction to Eulerian Grid
    call apply_full_correction(lag_grid,eul_grid,self%Cxcoef_local,self%Cxcoefo, &
        iside,steady_state_requested,normalization%dt_effective)
    
  end subroutine prepare_solver_step

  ! Set the previous step solution (frozen during GMRES iterations)
  ! This should be called ONCE at the start of each time step
  subroutine set_previous_step_solution(self, u_prev_in, dt_in)
    class(LinearSolver), intent(inout) :: self
    real(dp), intent(in) :: u_prev_in(:,:)
    real(dp), intent(in) :: dt_in
    
    ! Store frozen values - these won't change during GMRES iterations
    self%u_prev(1:nx,1:ny) = u_prev_in(1:nx,1:ny)
    self%dt_step = dt_in
    
  end subroutine set_previous_step_solution

  pure elemental logical function is_finite_backend_scalar(value) result(finite)
    real(dp),intent(in)::value
    integer(int64)::bits
    if(storage_size(value)==64 .and. radix(value)==2 .and. digits(value)==53 .and. &
        minexponent(value)==-1021 .and. maxexponent(value)==1024)then
      bits=transfer(value,bits)
      finite=ibits(bits,52,11)/=int(z'7ff',int64)
    else
      finite=ieee_is_finite(value)
    endif
  end function is_finite_backend_scalar

  pure elemental logical function safe_to_scale_backend(value,scale) result(safe)
    real(dp),intent(in)::value,scale
    integer::product_exponent
    integer(int64)::value_bits
    safe=.false.
    if(.not.is_finite_backend_scalar(value) .or. &
        .not.is_finite_backend_scalar(scale))return
    if(scale<=zero)return
    value_bits=transfer(value,value_bits)
    if(iand(value_bits,int(z'7fffffffffffffff',int64))==0_int64)then
      safe=.true.
      return
    endif
    product_exponent=exponent(abs(value))+exponent(scale)
    safe=product_exponent>minexponent(value) .and. &
        product_exponent<maxexponent(value)
  end function safe_to_scale_backend

  pure elemental logical function safe_to_divide_backend(value,divisor) result(safe)
    real(dp),intent(in)::value,divisor
    integer::quotient_exponent
    integer(int64)::value_bits
    safe=.false.
    if(.not.is_finite_backend_scalar(value) .or. &
        .not.is_finite_backend_scalar(divisor))return
    if(divisor<=zero)return
    value_bits=transfer(value,value_bits)
    if(iand(value_bits,int(z'7fffffffffffffff',int64))==0_int64)then
      safe=.true.
      return
    endif
    quotient_exponent=exponent(abs(value))-exponent(divisor)
    safe=quotient_exponent<maxexponent(value)
  end function safe_to_divide_backend

  pure elemental logical function safe_to_add_backend(left,right) result(safe)
    real(dp),intent(in)::left,right
    safe=.false.
    if(.not.is_finite_backend_scalar(left) .or. &
        .not.is_finite_backend_scalar(right))return
    if(left>zero .and. right>zero)then
      if(left>huge(one)-right)return
    elseif(left<zero .and. right<zero)then
      if(left<(-huge(one))-right)return
    endif
    safe=.true.
  end function safe_to_add_backend

  ! Active pump coefficient function
  ! Computes the spatially-varying jmp coefficient at arc-length s and time t
  pure function jmp_func(s, t) result(val)
    real(dp), intent(in) :: s, t
    real(dp) :: val

    ! PNAS SI Eq. (III.17), specialized to equal front/rear amplitudes and
    ! angular widths h_pl=t_pl=0.21*pi.
    val=0.0_dp
    if(t<=dualchem_pump_start_time+64.0_dp*epsilon(1.0_dp)* &
        max(1.0_dp,abs(dualchem_pump_start_time)))return
    val=pnas_gaussian_pump(s,kw)
  end function jmp_func

end module linear_solver_mod
