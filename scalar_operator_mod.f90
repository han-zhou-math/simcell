! Stage 10 typed contract for one semi-periodic scalar solve.
!
! Provenance: this narrow data contract is adapted from the independently
! developed dualchem scalar operator at commit 3b1407c.  Unlike the more
! general all-wall contract, this stage fixes the outer boundary topology:
! x is periodic and y has homogeneous total no-flux walls.  Four flux arrays
! remain in the typed layout for clarity, but all must be exact-zero sentinels
! because the inherited Stage 06 C++ entry point ignores wall-data arrays.
module scalar_operator_mod
  use ieee_arithmetic, only: ieee_is_finite
  use iso_fortran_env, only: int64
  use parameters, only: dp
  implicit none
  private

  integer, parameter, public :: SCALAR_STATUS_OK = 0
  integer, parameter, public :: SCALAR_STATUS_INVALID = -1

  type, public :: scalar_operator_t
    real(dp) :: diffusion = 1.0_dp
    real(dp) :: reaction = 0.0_dp
    real(dp), allocatable :: vx_face(:,:)
    real(dp), allocatable :: vy_face(:,:)
    real(dp), allocatable :: velocity_marker(:,:)
    real(dp), allocatable :: flux_w(:)
    real(dp), allocatable :: flux_e(:)
    real(dp), allocatable :: flux_s(:)
    real(dp), allocatable :: flux_n(:)
  end type scalar_operator_t

  type, public :: scalar_rhs_t
    real(dp), allocatable :: volume(:,:)
    real(dp), allocatable :: source_jump(:)
  end type scalar_rhs_t

  public :: validate_scalar_operator
  public :: validate_scalar_rhs
  public :: compute_scalar_lambda
  public :: normalize_scalar_system
  public :: is_finite_scalar

contains

  pure integer function validate_scalar_operator(operator, nx, ny, nmarkers) result(status)
    type(scalar_operator_t), intent(in) :: operator
    integer, intent(in) :: nx, ny, nmarkers

    status = SCALAR_STATUS_INVALID

    if (.not. valid_grid_extent(nx)) return
    if (.not. valid_grid_extent(ny)) return
    if (nmarkers < 0) return
    if (.not. is_finite_scalar(operator%diffusion)) return
    if (operator%diffusion <= 0.0_dp) return
    if (.not. is_finite_scalar(operator%reaction)) return
    if (operator%reaction < 0.0_dp) return

    if (.not. allocated(operator%vx_face)) return
    if (size(operator%vx_face, 1) /= nx + 1 .or. &
        size(operator%vx_face, 2) /= ny) return
    if (.not. allocated(operator%vy_face)) return
    if (size(operator%vy_face, 1) /= nx .or. &
        size(operator%vy_face, 2) /= ny + 1) return
    if (.not. allocated(operator%velocity_marker)) return
    if (size(operator%velocity_marker, 1) /= nmarkers .or. &
        size(operator%velocity_marker, 2) /= 2) return

    if (.not. allocated(operator%flux_w)) return
    if (size(operator%flux_w) /= ny) return
    if (.not. allocated(operator%flux_e)) return
    if (size(operator%flux_e) /= ny) return
    if (.not. allocated(operator%flux_s)) return
    if (size(operator%flux_s) /= nx) return
    if (.not. allocated(operator%flux_n)) return
    if (size(operator%flux_n) /= nx) return

    if (.not. all(is_finite_scalar(operator%vx_face))) return
    if (.not. all(is_finite_scalar(operator%vy_face))) return
    if (.not. all(is_finite_scalar(operator%velocity_marker))) return
    if (.not. all(is_finite_scalar(operator%flux_w))) return
    if (.not. all(is_finite_scalar(operator%flux_e))) return
    if (.not. all(is_finite_scalar(operator%flux_s))) return
    if (.not. all(is_finite_scalar(operator%flux_n))) return

    ! There is no west/east boundary in the semi-periodic operator.  Requiring
    ! exact zero catches an accidental all-wall call instead of silently
    ! discarding data at the periodic seam.
    if (any(operator%flux_w /= 0.0_dp)) return
    if (any(operator%flux_e /= 0.0_dp)) return
    if (any(operator%flux_s /= 0.0_dp)) return
    if (any(operator%flux_n /= 0.0_dp)) return

    status = SCALAR_STATUS_OK
  end function validate_scalar_operator

  pure integer function validate_scalar_rhs(rhs, nx, ny, nmarkers) result(status)
    type(scalar_rhs_t), intent(in) :: rhs
    integer, intent(in) :: nx, ny, nmarkers

    status = SCALAR_STATUS_INVALID
    if (.not. valid_grid_extent(nx)) return
    if (.not. valid_grid_extent(ny)) return
    if (nmarkers < 0) return
    if (.not. allocated(rhs%volume)) return
    if (size(rhs%volume, 1) /= nx .or. size(rhs%volume, 2) /= ny) return
    if (.not. allocated(rhs%source_jump)) return
    if (size(rhs%source_jump) /= nmarkers) return
    if (.not. all(is_finite_scalar(rhs%volume))) return
    if (.not. all(is_finite_scalar(rhs%source_jump))) return
    status = SCALAR_STATUS_OK
  end function validate_scalar_rhs

  pure subroutine compute_scalar_lambda(reaction, transient, dt, lambda, status)
    real(dp), intent(in) :: reaction, dt
    logical, intent(in) :: transient
    real(dp), intent(out) :: lambda
    integer, intent(out) :: status
    real(dp) :: inverse_dt, reciprocal_tiny, scaled_dt

    lambda = 0.0_dp
    status = SCALAR_STATUS_INVALID
    if (.not. is_finite_scalar(reaction)) return
    if (reaction < 0.0_dp) return

    if (.not. transient) then
      lambda = reaction
      status = SCALAR_STATUS_OK
      return
    endif

    if (.not. is_finite_scalar(dt)) return
    if (dt <= 0.0_dp) return

    ! Form 1/dt without trapping when dt is subnormal.  A nonrepresentable
    ! reciprocal is reported by status, not by an IEEE exception.
    if (dt < tiny(1.0_dp)) then
      reciprocal_tiny = 1.0_dp/tiny(1.0_dp)
      scaled_dt = dt/tiny(1.0_dp)
      if (scaled_dt < reciprocal_tiny/huge(1.0_dp)) return
      inverse_dt = reciprocal_tiny/scaled_dt
    else
      inverse_dt = 1.0_dp/dt
    endif
    if (.not. is_finite_scalar(inverse_dt)) return
    if (reaction > huge(1.0_dp) - inverse_dt) return

    lambda = reaction + inverse_dt
    if (.not. is_finite_scalar(lambda)) then
      lambda = 0.0_dp
      return
    endif
    status = SCALAR_STATUS_OK
  end subroutine compute_scalar_lambda

  pure subroutine normalize_scalar_system(operator, rhs, lambda_physical, &
       normalized_operator, normalized_rhs, lambda_normalized, status)
    type(scalar_operator_t), intent(in) :: operator
    type(scalar_rhs_t), intent(in) :: rhs
    real(dp), intent(in) :: lambda_physical
    type(scalar_operator_t), intent(out) :: normalized_operator
    type(scalar_rhs_t), intent(out) :: normalized_rhs
    real(dp), intent(out) :: lambda_normalized
    integer, intent(out) :: status
    real(dp) :: diffusion

    status = SCALAR_STATUS_INVALID
    lambda_normalized = 0.0_dp
    diffusion = operator%diffusion

    ! The caller validates shapes separately because nx, ny, and marker count
    ! belong to the frozen geometry.  This routine still fails safely if it is
    ! called with incomplete or nonfinite data.
    if (.not. is_finite_scalar(diffusion)) return
    if (diffusion <= 0.0_dp) return
    if (.not. is_finite_scalar(lambda_physical)) return
    if (lambda_physical < 0.0_dp) return
    if (.not. allocated(operator%vx_face)) return
    if (.not. allocated(operator%vy_face)) return
    if (.not. allocated(operator%velocity_marker)) return
    if (.not. allocated(operator%flux_w)) return
    if (.not. allocated(operator%flux_e)) return
    if (.not. allocated(operator%flux_s)) return
    if (.not. allocated(operator%flux_n)) return
    if (.not. allocated(rhs%volume)) return
    if (.not. allocated(rhs%source_jump)) return

    if (.not. all(safe_to_divide(operator%vx_face, diffusion))) return
    if (.not. all(safe_to_divide(operator%vy_face, diffusion))) return
    if (.not. all(safe_to_divide(operator%velocity_marker, diffusion))) return
    if (.not. all(safe_to_divide(operator%flux_w, diffusion))) return
    if (.not. all(safe_to_divide(operator%flux_e, diffusion))) return
    if (.not. all(safe_to_divide(operator%flux_s, diffusion))) return
    if (.not. all(safe_to_divide(operator%flux_n, diffusion))) return
    if (.not. all(safe_to_divide(rhs%volume, diffusion))) return
    if (.not. all(safe_to_divide(rhs%source_jump, diffusion))) return
    if (.not. safe_to_divide(operator%reaction, diffusion)) return
    if (.not. safe_to_divide(lambda_physical, diffusion)) return

    ! Intrinsic assignment performs a deep copy of the allocatable members.
    ! Only the disposable normalized objects are modified; the physical input
    ! contract remains immutable for residual callbacks and audit snapshots.
    normalized_operator = operator
    normalized_rhs = rhs
    normalized_operator%diffusion = 1.0_dp
    normalized_operator%reaction = operator%reaction/diffusion
    normalized_operator%vx_face = operator%vx_face/diffusion
    normalized_operator%vy_face = operator%vy_face/diffusion
    normalized_operator%velocity_marker = operator%velocity_marker/diffusion
    normalized_operator%flux_w = operator%flux_w/diffusion
    normalized_operator%flux_e = operator%flux_e/diffusion
    normalized_operator%flux_s = operator%flux_s/diffusion
    normalized_operator%flux_n = operator%flux_n/diffusion
    normalized_rhs%volume = rhs%volume/diffusion
    normalized_rhs%source_jump = rhs%source_jump/diffusion
    lambda_normalized = lambda_physical/diffusion

    if (.not. is_finite_scalar(lambda_normalized)) return
    status = SCALAR_STATUS_OK
  end subroutine normalize_scalar_system

  pure elemental logical function is_finite_scalar(value) result(is_finite)
    real(dp), intent(in) :: value
    integer(int64) :: bits

    ! ieee_is_finite can itself signal on a signaling NaN under the project's
    ! trap flags.  Inspect binary64 exponent bits on the supported platform.
    if (storage_size(value) == 64 .and. radix(value) == 2 .and. &
        digits(value) == 53 .and. minexponent(value) == -1021 .and. &
        maxexponent(value) == 1024) then
      bits = transfer(value, bits)
      is_finite = ibits(bits, 52, 11) /= int(z'7ff', int64)
    else
      is_finite = ieee_is_finite(value)
    endif
  end function is_finite_scalar

  pure logical function valid_grid_extent(extent) result(is_valid)
    integer, intent(in) :: extent
    is_valid = .false.
    if (extent < 2) return
    is_valid = iand(extent, extent - 1) == 0
  end function valid_grid_extent

  pure elemental logical function safe_to_divide(value, divisor) result(is_safe)
    real(dp), intent(in) :: value, divisor
    integer :: quotient_exponent

    is_safe = .false.
    if (.not. is_finite_scalar(value)) return
    if (.not. is_finite_scalar(divisor)) return
    if (divisor <= 0.0_dp) return
    if (value == 0.0_dp) then
      is_safe = .true.
      return
    endif

    ! This conservative exponent test is deliberately performed before the
    ! division so an extreme finite input cannot overflow under FPE traps.
    ! Equality is rejected because the significands decide the boundary case.
    quotient_exponent = exponent(abs(value)) - exponent(divisor)
    is_safe = quotient_exponent < maxexponent(value)
  end function safe_to_divide

end module scalar_operator_mod
