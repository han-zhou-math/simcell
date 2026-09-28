module grid_types
  use parameters, only: dp, nx, ny, xmin, xmax, ymin, ymax, hg, npts
  implicit none
  private
  public :: EulerianGrid, LagrangianGrid
  public :: init_eulerian, clean_eulerian, init_lagrangian, clean_lagrangian

  !-----------------------------------------------------------------------------
  ! Eulerian Grid Type
  !-----------------------------------------------------------------------------
  type :: EulerianGrid
    ! Grid dimensions and physical extent
    integer :: nx_grid, ny_grid
    real(dp) :: x_min, x_max, y_min, y_max
    real(dp) :: dx, dy

    ! Field variables
    real(dp), allocatable :: u(:,:), v(:,:), p(:,:)
    
    ! Geometry tags and distance maps
    ! id: 0=fluid, 1=solid/interface
    ! idf: face tags
    integer, allocatable :: id(:,:), idf(:,:)
    
    ! Distance maps and closest IB point indices
    ! dmapc: distance map current, dmapo: distance map old
    real(dp), allocatable :: dmapc(:,:), dmapo(:,:) 
    
    ! kaic: index of closest IB point (current)
    ! kaio: index of closest IB point (old)
    integer, allocatable :: kaic(:,:), kaio(:,:)
    
    ! Correction function (for C++ solver)
    real(dp), allocatable :: crc(:,:)
    
    ! Auxiliary arrays for tracking status
    integer, allocatable :: chkero(:,:), chkerc(:,:)
    integer, allocatable :: oid(:,:), qid(:,:)
    integer, allocatable :: idn(:,:) ! Freshly cleared tags

  contains
    procedure :: init => init_eulerian
    procedure :: clean => clean_eulerian
  end type EulerianGrid

  !-----------------------------------------------------------------------------
  ! Lagrangian Grid Type
  !-----------------------------------------------------------------------------
  type :: LagrangianGrid
    ! Number of Lagrangian points (e.g., npts)
    integer :: npts
    
    ! Coordinates
    real(dp), allocatable :: x(:), y(:)
    real(dp), allocatable :: x_old(:), y_old(:)
    
    ! Geometry information
    real(dp), allocatable :: normal(:,:)  ! (npts, 2)
    real(dp), allocatable :: tangent(:,:) ! (npts, 2)
    real(dp), allocatable :: mk(:,:)      ! (7, npts) - Polynomial coeffs
    
    ! Geometric constraints / Communication links
    ! Stores indices of the Eulerian grid cell containing the IB point
    ! Corresponds to `plinkij` in current code
    ! (npts, 4): [i, j, i_corner, j_corner]
    ! Stores the indices of the Eulerian grid cell (bottom-left corner)
    ! associated with each IB point.
    integer, allocatable :: plinkij(:,:) 
    
    ! Solver variables (defined at IB points)
    ! For GMRES, this holds the vector 'x' or 'b'
    ! Scalar variable per point
    real(dp), allocatable :: vars(:) 
    
    ! BC information and exact solution
    real(dp), allocatable :: rhs(:)  ! For storing get_u_BC
    real(dp), allocatable :: phi(:)  ! For storing get_u_exact 

  contains
    procedure :: init => init_lagrangian
    procedure :: clean => clean_lagrangian
    procedure :: advance => advance_lagrangian
    procedure :: check_bounds => check_lagrangian_bounds
    
    ! Communication procedures for Matrix-Free Product
    procedure :: map_to_eulerian 
    procedure :: map_from_eulerian
  end type LagrangianGrid

contains

  !-----------------------------------------------------------------------------
  ! Eulerian Grid Methods
  !-----------------------------------------------------------------------------
  subroutine init_eulerian(self, nx_in, ny_in, xmin_in, xmax_in, ymin_in, ymax_in)
    class(EulerianGrid), intent(inout) :: self
    integer, intent(in), optional :: nx_in, ny_in
    real(dp), intent(in), optional :: xmin_in, xmax_in, ymin_in, ymax_in
    
    if (present(nx_in)) then
        self%nx_grid = nx_in
    else
        self%nx_grid = nx
    endif
    
    if (present(ny_in)) then
        self%ny_grid = ny_in
    else
        self%ny_grid = ny
    endif
    
    if (present(xmin_in)) then
        self%x_min = xmin_in
    else
        self%x_min = xmin
    endif
    
    if (present(xmax_in)) then
        self%x_max = xmax_in
    else
        self%x_max = xmax
    endif
    
    if (present(ymin_in)) then
        self%y_min = ymin_in
    else
        self%y_min = ymin
    endif
    
    if (present(ymax_in)) then
        self%y_max = ymax_in
    else
        self%y_max = ymax
    endif

    self%dx = (self%x_max - self%x_min) / real(self%nx_grid, dp)
    self%dy = (self%y_max - self%y_min) / real(self%ny_grid, dp)
    
    ! Check consistency with global hg if using globals? 
    ! For now, just calculate dx/dy from bounds.
    ! But wait, hg is used elsewhere. If we pass custom bounds, hg might be wrong?
    ! The user passed 0.0, 1.0, 0.0, 1.0 with nx, ny.
    ! If nx=128, dx = 1/128.
    ! Global hg might be different.
    ! However, for the test, it should be fine.

    ! Allocate arrays with ghost cells if needed, matching original code (-1:nx+1, -1:ny+1)
    ! Use local dimensions
    allocate(self%u(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%v(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%p(-1:self%nx_grid+1, -1:self%ny_grid+1))
    
    allocate(self%id(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%idf(-1:self%nx_grid+1, -1:self%ny_grid+1))
    
    allocate(self%dmapc(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%dmapo(-1:self%nx_grid+1, -1:self%ny_grid+1))
    
    allocate(self%kaic(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%kaio(-1:self%nx_grid+1, -1:self%ny_grid+1))
    
    allocate(self%crc(self%nx_grid, self%ny_grid))
    
    allocate(self%chkero(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%chkerc(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%oid(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%qid(-1:self%nx_grid+1, -1:self%ny_grid+1))
    allocate(self%idn(-1:self%nx_grid+1, -1:self%ny_grid+1))

    ! Initialize to zero
    self%u = 0.0_dp
    self%v = 0.0_dp
    self%p = 0.0_dp
    self%id = 0
    self%idf = 0
    self%dmapc = 0.0_dp
    self%dmapo = 0.0_dp
    self%kaic = 0
    self%kaio = 0
    self%crc = 0.0_dp
    self%chkero = 0
    self%chkerc = 0
    self%oid = 0
    self%qid = 0
    self%idn = 0

  end subroutine init_eulerian

  subroutine clean_eulerian(self)
    class(EulerianGrid), intent(inout) :: self
    if (allocated(self%u)) deallocate(self%u)
    if (allocated(self%v)) deallocate(self%v)
    if (allocated(self%p)) deallocate(self%p)
    if (allocated(self%id)) deallocate(self%id)
    if (allocated(self%idf)) deallocate(self%idf)
    if (allocated(self%dmapc)) deallocate(self%dmapc)
    if (allocated(self%dmapo)) deallocate(self%dmapo)
    if (allocated(self%kaic)) deallocate(self%kaic)
    if (allocated(self%kaio)) deallocate(self%kaio)
    if (allocated(self%crc)) deallocate(self%crc)
    if (allocated(self%chkero)) deallocate(self%chkero)
    if (allocated(self%chkerc)) deallocate(self%chkerc)
    if (allocated(self%oid)) deallocate(self%oid)
    if (allocated(self%qid)) deallocate(self%qid)
    if (allocated(self%idn)) deallocate(self%idn)
  end subroutine clean_eulerian

  !-----------------------------------------------------------------------------
  ! Lagrangian Grid Methods
  !-----------------------------------------------------------------------------
  subroutine init_lagrangian(self, num_pts)
    class(LagrangianGrid), intent(inout) :: self
    integer, intent(in) :: num_pts
    
    self%npts = num_pts
    
    allocate(self%x(num_pts))
    allocate(self%y(num_pts))
    allocate(self%x_old(num_pts))
    allocate(self%y_old(num_pts))
    
    allocate(self%normal(num_pts, 2))
    allocate(self%tangent(num_pts, 2))
    allocate(self%mk(7, num_pts))
    
    allocate(self%plinkij(num_pts, 4))
    
    ! Initialize vars with a default size, can be reallocated if needed
    ! Initialize vars with a default size, can be reallocated if needed
    allocate(self%vars(num_pts)) 
    allocate(self%rhs(num_pts))
    allocate(self%phi(num_pts)) 

    self%x = 0.0_dp
    self%y = 0.0_dp
    self%x_old = 0.0_dp
    self%y_old = 0.0_dp
    self%normal = 0.0_dp
    self%tangent = 0.0_dp
    self%mk = 0.0_dp
!    self%plinkij = 0
    self%plinkij = 0
    self%vars = 0.0_dp
    self%rhs = 0.0_dp
    self%phi = 0.0_dp

  end subroutine init_lagrangian

  subroutine clean_lagrangian(self)
    class(LagrangianGrid), intent(inout) :: self
    if (allocated(self%x)) deallocate(self%x)
    if (allocated(self%y)) deallocate(self%y)
    if (allocated(self%x_old)) deallocate(self%x_old)
    if (allocated(self%y_old)) deallocate(self%y_old)
    if (allocated(self%normal)) deallocate(self%normal)
    if (allocated(self%tangent)) deallocate(self%tangent)
    if (allocated(self%mk)) deallocate(self%mk)
    if (allocated(self%plinkij)) deallocate(self%plinkij)
    if (allocated(self%vars)) deallocate(self%vars)
    if (allocated(self%rhs)) deallocate(self%rhs)
    if (allocated(self%phi)) deallocate(self%phi)
  end subroutine clean_lagrangian

  subroutine advance_lagrangian(self, dt, vel_ib)
    class(LagrangianGrid), intent(inout) :: self
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: vel_ib(:,:)  ! (npts, 2) - velocity at IB points
    
    integer :: k
    
    ! Save current positions as old
    self%x_old = self%x
    self%y_old = self%y
    
    ! Advance positions using velocity
    do k = 1, self%npts
      self%x(k) = self%x(k) + dt * vel_ib(k, 1)
      self%y(k) = self%y(k) + dt * vel_ib(k, 2)
    end do
    
  end subroutine advance_lagrangian

  subroutine check_lagrangian_bounds(self, x_lo, x_hi, y_lo, y_hi, istep)
    ! Check that all Lagrangian points are within the specified bounds
    ! Stops simulation with error message if any point is out of bounds
    class(LagrangianGrid), intent(in) :: self
    real(dp), intent(in) :: x_lo, x_hi, y_lo, y_hi
    integer, intent(in), optional :: istep
    
    integer :: k, step_num
    logical :: out_of_bounds
    
    if (present(istep)) then
        step_num = istep
    else
        step_num = -1
    endif
    
    out_of_bounds = .false.
    
    do k = 1, self%npts
        if (self%x(k) < x_lo .or. self%x(k) > x_hi .or. &
            self%y(k) < y_lo .or. self%y(k) > y_hi) then
            print *, "***********************************************************"
            print *, "ERROR: Interface point out of bounds"
            print *, "  Point index:", k
            if (step_num >= 0) print *, "  Time step:", step_num
            print *, "  Position: (", self%x(k), ",", self%y(k), ")"
            print *, "  Valid range: x in [", x_lo, ",", x_hi, "]"
            print *, "              y in [", y_lo, ",", y_hi, "]"
            print *, "***********************************************************"
            out_of_bounds = .true.
        endif
    enddo
    
    if (out_of_bounds) then
        print *, "Stopping simulation: interface too close to domain boundary"
        stop
    endif
    
  end subroutine check_lagrangian_bounds

  ! Standard 4-point Delta Function
  function delta_func(r) result(val)
    real(dp), intent(in) :: r
    real(dp) :: val, abs_r
    
    abs_r = abs(r)
    if (abs_r <= 1.0_dp) then
        val = 0.125_dp * (3.0_dp - 2.0_dp*abs_r + sqrt(1.0_dp + 4.0_dp*abs_r - 4.0_dp*abs_r*abs_r))
    elseif (abs_r <= 2.0_dp) then
        val = 0.125_dp * (5.0_dp - 2.0_dp*abs_r - sqrt(-7.0_dp + 12.0_dp*abs_r - 4.0_dp*abs_r*abs_r))
    else
        val = 0.0_dp
    endif
  end function delta_func

  subroutine map_to_eulerian(self, eul_grid)
    class(LagrangianGrid), intent(in) :: self
    type(EulerianGrid), intent(inout) :: eul_grid
    
    integer :: k, i, j, i_min, i_max, j_min, j_max
    real(dp) :: x, y, dx, dy, val, factor
    real(dp) :: w_x, w_y
    
    eul_grid%u = 0.0_dp ! Clear target field (assuming we map to u)
    
    dx = eul_grid%dx
    dy = eul_grid%dy
    
    ! Scaling factor for spreading force: usually 1/h^2?
    ! Or just sum(w * F).
    ! If p is force density on boundary (force per unit length), then F_grid = sum(p * delta * ds).
    ! Here self%vars is p.
    ! We assume p includes the ds factor or we need to multiply by ds.
    ! Let's assume p is the force at the point (Lagrangian force).
    ! Then F_grid(x) = sum_k p_k * delta(x - x_k).
    ! Note: delta has units 1/length^2.
    ! So F_grid has units Force / length^2 = Force density.
    ! If p_k is Force * length (integrated force?), no p_k is Force.
    ! Wait, usually F(x) = int F(s) delta(x - X(s)) ds.
    ! Discretized: F_ij = sum_k F_k * delta_h(x_ij - X_k) * ds_k.
    ! self%vars usually stores F_k * ds_k or just F_k?
    ! In standard IB, we store Force F_k. And multiply by ds_k during spreading.
    ! Here, let's assume self%vars is just the coefficient to be spread.
    ! We need ds_k.
    ! For uniform circle, ds = 2*pi*r / npts.
    ! Let's compute ds on the fly or assume uniform.
    ! For now, let's assume self%vars includes the weight (p * ds).
    ! Or just spread self%vars directly.
    ! The linear solver solves A p = b.
    ! If A includes spreading and interpolation, and we want symmetric A,
    ! S (spread) and S^T (interp) should be adjoint.
    ! S: L -> E. (Spread)
    ! S^T: E -> L. (Interpolate)
    ! Interpolate: U_k = sum_ij u_ij * delta_h(x_ij - X_k) * hx * hy.
    ! Spread: f_ij = sum_k F_k * delta_h(x_ij - X_k).
    ! Are these adjoint?
    ! <S F, u>_E = sum_ij (sum_k F_k delta) u_ij hx hy = sum_k F_k (sum_ij u_ij delta hx hy) = <F, S^T u>_L.
    ! Yes, if inner product on E is weighted by hx*hy, and on L is just sum (unweighted).
    ! So:
    ! Spread: f_ij = sum_k F_k * delta(x_ij - X_k)
    ! Interp: U_k = sum_ij u_ij * delta(x_ij - X_k) * hx * hy
    
    ! Let's implement this.
    
    do k = 1, self%npts
        x = self%x(k)
        y = self%y(k)
        
        ! Check if point is within domain bounds
        if (x < eul_grid%x_min .or. x > eul_grid%x_max .or. &
            y < eul_grid%y_min .or. y > eul_grid%y_max) then
            print *, '***********************************************************'
            print *, 'ERROR: Lagrangian point outside domain in map_to_eulerian'
            print *, '  Point index k=', k
            print *, '  Position (x,y)=', x, y
            print *, '  Domain bounds: x=[', eul_grid%x_min, ',', eul_grid%x_max, ']'
            print *, '                 y=[', eul_grid%y_min, ',', eul_grid%y_max, ']'
            print *, '***********************************************************'
            stop 'Lagrangian point outside domain'
        endif
        
        val = self%vars(k)
        
        ! Find grid range (support is +/- 2h)
        i_min = int((x - eul_grid%x_min) / dx - 2.0_dp)
        i_max = int((x - eul_grid%x_min) / dx + 3.0_dp)
        j_min = int((y - eul_grid%y_min) / dy - 2.0_dp)
        j_max = int((y - eul_grid%y_min) / dy + 3.0_dp)
        
        do j = j_min, j_max
            do i = i_min, i_max
                ! Check bounds (periodic or clamp? Assume clamp or ignore if out)
                if (i >= 1 .and. i <= eul_grid%nx_grid .and. j >= 1 .and. j <= eul_grid%ny_grid) then
                    w_x = delta_func((x - (eul_grid%x_min + (real(i,dp)-0.5_dp)*dx)) / dx)
                    w_y = delta_func((y - (eul_grid%y_min + (real(j,dp)-0.5_dp)*dy)) / dy)
                    
                    eul_grid%u(i, j) = eul_grid%u(i, j) + val * w_x * w_y / (dx * dy) 
                    ! Note: delta_func returns phi(r). delta = phi(x/h)/h.
                    ! So delta_2d = phi(rx)/dx * phi(ry)/dy.
                    ! So we divide by dx*dy.
                endif
            enddo
        enddo
    enddo
  end subroutine map_to_eulerian

  subroutine map_from_eulerian(self, eul_grid)
    class(LagrangianGrid), intent(inout) :: self
    type(EulerianGrid), intent(in) :: eul_grid
    
    integer :: k, i, j, i_min, i_max, j_min, j_max
    real(dp) :: x, y, dx, dy, val
    real(dp) :: w_x, w_y
    
    dx = eul_grid%dx
    dy = eul_grid%dy
    
    do k = 1, self%npts
        x = self%x(k)
        y = self%y(k)
        
        ! Check if point is within domain bounds
        if (x < eul_grid%x_min .or. x > eul_grid%x_max .or. &
            y < eul_grid%y_min .or. y > eul_grid%y_max) then
            print *, '***********************************************************'
            print *, 'ERROR: Lagrangian point outside domain in map_from_eulerian'
            print *, '  Point index k=', k
            print *, '  Position (x,y)=', x, y
            print *, '  Domain bounds: x=[', eul_grid%x_min, ',', eul_grid%x_max, ']'
            print *, '                 y=[', eul_grid%y_min, ',', eul_grid%y_max, ']'
            print *, '***********************************************************'
            stop 'Lagrangian point outside domain'
        endif
        
        val = 0.0_dp
        
        i_min = int((x - eul_grid%x_min) / dx - 2.0_dp)
        i_max = int((x - eul_grid%x_min) / dx + 3.0_dp)
        j_min = int((y - eul_grid%y_min) / dy - 2.0_dp)
        j_max = int((y - eul_grid%y_min) / dy + 3.0_dp)
        
        do j = j_min, j_max
            do i = i_min, i_max
                if (i >= 1 .and. i <= eul_grid%nx_grid .and. j >= 1 .and. j <= eul_grid%ny_grid) then
                    w_x = delta_func((x - (eul_grid%x_min + (real(i,dp)-0.5_dp)*dx)) / dx)
                    w_y = delta_func((y - (eul_grid%y_min + (real(j,dp)-0.5_dp)*dy)) / dy)
                    
                    ! Interpolation: sum u_ij * delta * hx * hy
                    ! delta = phi/h.
                    ! val += u_ij * (phi_x/dx * phi_y/dy) * dx * dy
                    ! val += u_ij * phi_x * phi_y
                    
                    val = val + eul_grid%u(i, j) * w_x * w_y
                endif
            enddo
        enddo
        self%vars(k) = val
    enddo
  end subroutine map_from_eulerian

end module grid_types
