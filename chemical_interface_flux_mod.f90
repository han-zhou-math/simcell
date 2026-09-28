! SPDX-License-Identifier: BSD-3-Clause
!
! Pure moving-interface chemical-flux algebra for Stage 06 (Stage 6A port,
! hardened in the free-actin style).
!
! At a marker the one-sided Robin law is written in physical units as
!
!     -D * dC/dn = u*C + o*C_opp + s,
!
! where C is the trace on the side being solved and C_opp the opposite-side
! trace (the cross-membrane coupling that free actin does not have).  With
! the relative normal velocity
!
!     r = (v_marker - V_Gamma) . n,     V_Gamma = (X^{n+1}-X^n)/dt,
!
! and the split pump coefficients
!
!     a = kc + max(p,0),   b = -kc + min(p,0),
!
! the four-tuple is
!
!     exterior:  u = (r-b)/D,   o = a/D,   s = (p-r)/D,   gradient = -1
!     interior:  u = (r-a)/D,   o = b/D,   s = (p-r)/D,   gradient = -1
!
! The diffusion coefficient is carried explicitly (the free-actin convention
! of dividing by D) instead of being asserted to one, and the coefficient
! kernel is pure elemental so each marker is independently testable.  The
! exterior/interior asymmetry in a and b and the opposite-trace coupling o
! are preserved exactly; at D=1 the coefficients reduce bit-for-bit to the
! validated Stage 6A form.
module chemical_interface_flux_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use parameters, only: dp, one, zero
  implicit none
  private

  integer, parameter, public :: CHEMICAL_INTERIOR = 1
  integer, parameter, public :: CHEMICAL_EXTERIOR = -1
  integer, parameter, public :: CHEMICAL_FLUX_OK = 0
  integer, parameter, public :: CHEMICAL_FLUX_INVALID = 1

  type, public :: chemical_normalization_t
    real(dp) :: diffusion = 0.0_dp
    real(dp) :: inverse_diffusion = 0.0_dp
    real(dp) :: dt_effective = 0.0_dp
  end type chemical_normalization_t

  public :: build_chemical_normalization
  public :: normalize_chemical_absolute_tolerance
  public :: build_relative_normal_velocity
  public :: build_chemical_robin

contains

  pure subroutine build_chemical_normalization(diffusion,dt,normalization,status)
    real(dp), intent(in) :: diffusion,dt
    type(chemical_normalization_t), intent(out) :: normalization
    integer, intent(out) :: status

    call reset_normalization(normalization)
    status=CHEMICAL_FLUX_INVALID
    if(.not.is_finite_scalar(diffusion) .or. .not.is_finite_scalar(dt))return
    if(diffusion<=zero .or. dt<=zero)return
    if(.not.safe_to_divide(one,diffusion))return
    if(.not.safe_to_multiply(diffusion,dt))return

    normalization%diffusion=diffusion
    normalization%inverse_diffusion=one/diffusion
    normalization%dt_effective=diffusion*dt
    if(.not.is_finite_scalar(normalization%inverse_diffusion) .or. &
        normalization%inverse_diffusion<=zero .or. &
        .not.is_finite_scalar(normalization%dt_effective) .or. &
        normalization%dt_effective<=zero)then
      call reset_normalization(normalization)
      return
    endif
    status=CHEMICAL_FLUX_OK
  end subroutine build_chemical_normalization

  pure subroutine normalize_chemical_absolute_tolerance(physical_atol, &
      diffusion,normalized_atol,status)
    real(dp),intent(in) :: physical_atol,diffusion
    real(dp),intent(out) :: normalized_atol
    integer,intent(out) :: status

    normalized_atol=zero
    status=CHEMICAL_FLUX_INVALID
    if(.not.is_finite_scalar(physical_atol) .or. physical_atol<zero)return
    if(.not.safe_to_divide(physical_atol,diffusion))return
    normalized_atol=physical_atol/diffusion
    if(.not.is_finite_scalar(normalized_atol) .or. normalized_atol<zero)then
      normalized_atol=zero
      return
    endif
    status=CHEMICAL_FLUX_OK
  end subroutine normalize_chemical_absolute_tolerance

  pure subroutine build_relative_normal_velocity(marker_velocity,interface_velocity, &
      normal,r,status)
    real(dp), intent(in) :: marker_velocity(:,:), interface_velocity(:,:), normal(:,:)
    real(dp), intent(out) :: r(:)
    integer, intent(out) :: status
    real(dp) :: normal_length(size(normal,1))
    real(dp) :: tolerance

    r = zero
    status = CHEMICAL_FLUX_INVALID
    if (size(marker_velocity,2) /= 2 .or. size(interface_velocity,2) /= 2 .or. &
        size(normal,2) /= 2) return
    if (size(marker_velocity,1) <= 0) return
    if (size(interface_velocity,1) /= size(marker_velocity,1) .or. &
        size(normal,1) /= size(marker_velocity,1) .or. &
        size(r) /= size(marker_velocity,1)) return
    if (.not. all(ieee_is_finite(marker_velocity)) .or. &
        .not. all(ieee_is_finite(interface_velocity)) .or. &
        .not. all(ieee_is_finite(normal))) return

    normal_length = sqrt(sum(normal**2,dim=2))
    tolerance = 256.0_dp*epsilon(one)
    if (maxval(abs(normal_length-one)) > tolerance) return

    r = sum((marker_velocity-interface_velocity)*normal,dim=2)
    if (.not. all(ieee_is_finite(r))) then
      r = zero
      return
    end if
    status = CHEMICAL_FLUX_OK
  end subroutine build_relative_normal_velocity

  pure elemental subroutine build_chemical_robin(side,r,kc,p,diffusion, &
      unknown_coef,gradient_coef,opposite_coef,shift_constant,status)
    integer, intent(in) :: side
    real(dp), intent(in) :: r, kc, p, diffusion
    real(dp), intent(out) :: unknown_coef, gradient_coef, opposite_coef, shift_constant
    integer, intent(out) :: status
    real(dp) :: a,b,b_magnitude,p_minus,p_plus
    real(dp) :: raw_unknown,raw_opposite,raw_shift

    call reset_robin(unknown_coef,gradient_coef,opposite_coef,shift_constant)
    status = CHEMICAL_FLUX_INVALID
    if (.not. is_finite_scalar(r) .or. .not. is_finite_scalar(kc) .or. &
        .not. is_finite_scalar(p) .or. .not. is_finite_scalar(diffusion)) return
    if (kc < zero .or. diffusion <= zero) return

    p_plus = max(p,zero)
    p_minus = min(p,zero)
    if(.not.safe_to_add(kc,p_plus))return
    a = kc+p_plus
    if(.not.safe_to_add(kc,-p_minus))return
    b_magnitude=kc-p_minus
    b=-b_magnitude
    select case(side)
    case(CHEMICAL_EXTERIOR)
      if(.not.safe_to_subtract(r,b))return
      raw_unknown=r-b
      raw_opposite=a
    case(CHEMICAL_INTERIOR)
      if(.not.safe_to_subtract(r,a))return
      raw_unknown=r-a
      raw_opposite=b
    case default
      return
    end select
    if(.not.safe_to_subtract(p,r))return
    raw_shift=p-r
    if(.not.safe_to_divide(raw_unknown,diffusion) .or. &
        .not.safe_to_divide(raw_opposite,diffusion) .or. &
        .not.safe_to_divide(raw_shift,diffusion))return
    unknown_coef=raw_unknown/diffusion
    gradient_coef = -one
    opposite_coef=raw_opposite/diffusion
    shift_constant=raw_shift/diffusion
    if(.not.is_finite_scalar(unknown_coef) .or. &
        .not.is_finite_scalar(gradient_coef) .or. &
        .not.is_finite_scalar(opposite_coef) .or. &
        .not.is_finite_scalar(shift_constant))then
      call reset_robin(unknown_coef,gradient_coef,opposite_coef,shift_constant)
      return
    endif
    status = CHEMICAL_FLUX_OK
  end subroutine build_chemical_robin

  pure subroutine reset_normalization(normalization)
    type(chemical_normalization_t),intent(out)::normalization
    normalization%diffusion=zero
    normalization%inverse_diffusion=zero
    normalization%dt_effective=zero
  end subroutine reset_normalization

  pure elemental subroutine reset_robin(unknown_coef,gradient_coef,opposite_coef,shift_constant)
    real(dp),intent(out)::unknown_coef,gradient_coef,opposite_coef,shift_constant
    unknown_coef=zero
    gradient_coef=zero
    opposite_coef=zero
    shift_constant=zero
  end subroutine reset_robin

  pure elemental logical function is_finite_scalar(value) result(is_finite)
    real(dp),intent(in)::value
    integer(int64)::bits
    if(storage_size(value)==64 .and. radix(value)==2 .and. digits(value)==53 .and. &
        minexponent(value)==-1021 .and. maxexponent(value)==1024)then
      bits=transfer(value,bits)
      is_finite=ibits(bits,52,11)/=int(z'7ff',int64)
    else
      is_finite=ieee_is_finite(value)
    endif
  end function is_finite_scalar

  pure elemental logical function safe_to_divide(value,divisor)
    real(dp), intent(in) :: value,divisor
    integer :: quotient_exponent
    integer(int64)::value_bits

    safe_to_divide=.false.
    if(.not.is_finite_scalar(value) .or. .not.is_finite_scalar(divisor))return
    if(divisor<=zero)return
    value_bits=transfer(value,value_bits)
    if(iand(value_bits,int(z'7fffffffffffffff',int64))==0_int64)then
      safe_to_divide=.true.
      return
    endif
    quotient_exponent=exponent(abs(value))-exponent(divisor)
    safe_to_divide=quotient_exponent<maxexponent(value)
  end function safe_to_divide

  pure elemental logical function safe_to_multiply(left,right)
    real(dp), intent(in) :: left,right
    integer :: product_exponent

    safe_to_multiply=.false.
    if(.not.is_finite_scalar(left) .or. .not.is_finite_scalar(right))return
    if(left<=zero .or. right<=zero)return
    product_exponent=exponent(left)+exponent(right)
    if(product_exponent<=minexponent(left) .or. &
        product_exponent>=maxexponent(left))return
    safe_to_multiply=.true.
  end function safe_to_multiply

  pure elemental logical function safe_to_add(left,right)
    real(dp),intent(in)::left,right
    safe_to_add=.false.
    if(.not.is_finite_scalar(left) .or. .not.is_finite_scalar(right))return
    if(left>zero .and. right>zero)then
      if(left>huge(one)-right)return
    elseif(left<zero .and. right<zero)then
      if(left < -huge(one)-right)return
    endif
    safe_to_add=.true.
  end function safe_to_add

  pure elemental logical function safe_to_subtract(left,right)
    real(dp),intent(in)::left,right
    safe_to_subtract=safe_to_add(left,-right)
  end function safe_to_subtract

end module chemical_interface_flux_mod
