! Stage 10 reduced two-actin coefficient and boundary algebra.
!
! Provenance: adapted narrowly from code/dualchem commit bac71cc.  That source
! supplied the overflow-safe evaluation of a=eta/(eta+eta_s) and
! D_n=k_sigma/(eta+eta_s).  Stage 10 adds explicit normalized Robin helpers so
! the network off-diagonal sign and the unscaled interface velocity are visible
! without invoking a Cartesian or interface solver.
module actin_model_mod
  use ieee_arithmetic, only: ieee_is_finite
  use iso_fortran_env, only: int64
  use parameters, only: dp
  implicit none
  private

  integer, parameter, public :: ACTIN_OK = 0
  integer, parameter, public :: ACTIN_INVALID_COEFFICIENTS = 1

  type, public :: actin_coefficients_t
    real(dp) :: network_transport_fraction = 0.0_dp
    real(dp) :: network_diffusivity = 0.0_dp
    real(dp) :: free_diffusivity = 0.0_dp
    real(dp) :: turnover_rate = 0.0_dp
    real(dp) :: membrane_rate = 0.0_dp
  end type actin_coefficients_t

  public :: initialize_actin_coefficients
  public :: network_boundary_residual, free_boundary_residual
  public :: network_robin_residual, free_robin_residual
  public :: network_robin_beta, free_robin_beta, membrane_coupling_q

contains

  pure subroutine initialize_actin_coefficients(eta, eta_s, k_sigma, d_c, &
       gamma, j_c, coeff, status)
    real(dp), intent(in) :: eta, eta_s, k_sigma, d_c, gamma, j_c
    type(actin_coefficients_t), intent(out) :: coeff
    integer, intent(out) :: status
    integer :: denominator_exponent, quotient_exponent
    real(dp) :: eta_fraction, eta_s_fraction, denominator_fraction
    real(dp) :: transport_ratio, quotient, quotient_fraction
    real(dp) :: network_transport_fraction, network_diffusivity

    coeff = actin_coefficients_t()
    status = ACTIN_INVALID_COEFFICIENTS

    ! The bit-level finite predicate below rejects a signaling NaN without
    ! evaluating a floating-point comparison under -ffpe-trap=invalid.
    if (.not. is_finite_coefficient(eta)) return
    if (.not. is_finite_coefficient(eta_s)) return
    if (.not. is_finite_coefficient(k_sigma)) return
    if (.not. is_finite_coefficient(d_c)) return
    if (.not. is_finite_coefficient(gamma)) return
    if (.not. is_finite_coefficient(j_c)) return

    if (eta <= 0.0_dp .or. eta_s <= 0.0_dp .or. k_sigma <= 0.0_dp) return
    if (d_c <= 0.0_dp .or. gamma < 0.0_dp .or. j_c < 0.0_dp) return

    ! Scale the denominator before addition, so eta+eta_s never overflows.
    denominator_exponent = max(exponent(eta), exponent(eta_s))
    eta_fraction = scale(fraction(eta), exponent(eta)-denominator_exponent)
    eta_s_fraction = scale(fraction(eta_s), &
         exponent(eta_s)-denominator_exponent)
    denominator_fraction = eta_fraction + eta_s_fraction

    ! Ordered ratio evaluation retains a representable minimum-subnormal a.
    if (eta < eta_s) then
      transport_ratio = eta/eta_s
      network_transport_fraction = transport_ratio/(1.0_dp+transport_ratio)
    else
      transport_ratio = eta_s/eta
      network_transport_fraction = 1.0_dp/(1.0_dp+transport_ratio)
    endif

    ! Reconstruct k_sigma/(eta+eta_s) from fractions and exponents.  This
    ! avoids both huge+huge overflow and a staged-underflow loss.
    quotient = fraction(k_sigma)/denominator_fraction
    quotient_fraction = fraction(quotient)
    quotient_exponent = exponent(quotient)+exponent(k_sigma) &
         -denominator_exponent
    if (quotient_exponent > maxexponent(1.0_dp)) return
    if (quotient_exponent == maxexponent(1.0_dp)) then
      if (quotient_fraction > fraction(huge(1.0_dp))) return
    endif

    network_diffusivity = scale(quotient_fraction, quotient_exponent)
    if (.not. is_finite_coefficient(network_diffusivity)) return
    if (network_diffusivity <= 0.0_dp) return

    coeff%network_transport_fraction = network_transport_fraction
    coeff%network_diffusivity = network_diffusivity
    coeff%free_diffusivity = d_c
    coeff%turnover_rate = gamma
    coeff%membrane_rate = j_c
    status = ACTIN_OK
  end subroutine initialize_actin_coefficients

  pure elemental function network_boundary_residual(coeff, theta_n, theta_c, &
       vc_dot_n, interface_velocity_dot_n, grad_theta_n_dot_n) result(residual)
    type(actin_coefficients_t), intent(in) :: coeff
    real(dp), intent(in) :: theta_n, theta_c, vc_dot_n
    real(dp), intent(in) :: interface_velocity_dot_n, grad_theta_n_dot_n
    real(dp) :: residual

    ! ((a*v_c-V_Gamma)*theta_n-D_n*grad(theta_n)).n=-j_c*theta_c
    residual = coeff%network_transport_fraction*theta_n*vc_dot_n &
         -coeff%network_diffusivity*grad_theta_n_dot_n &
         -theta_n*interface_velocity_dot_n+coeff%membrane_rate*theta_c
  end function network_boundary_residual

  pure elemental function free_boundary_residual(coeff, theta_c, vc_dot_n, &
       interface_velocity_dot_n, grad_theta_c_dot_n) result(residual)
    type(actin_coefficients_t), intent(in) :: coeff
    real(dp), intent(in) :: theta_c, vc_dot_n
    real(dp), intent(in) :: interface_velocity_dot_n, grad_theta_c_dot_n
    real(dp) :: residual

    ! ((v_c-V_Gamma)*theta_c-D_c*grad(theta_c)).n=j_c*theta_c
    residual = theta_c*(vc_dot_n-interface_velocity_dot_n) &
         -coeff%free_diffusivity*grad_theta_c_dot_n &
         -coeff%membrane_rate*theta_c
  end function free_boundary_residual

  pure elemental function network_robin_beta(coeff, vc_dot_n, &
       interface_velocity_dot_n) result(beta)
    type(actin_coefficients_t), intent(in) :: coeff
    real(dp), intent(in) :: vc_dot_n, interface_velocity_dot_n
    real(dp) :: beta

    ! Only v_c is scaled by a.  V_Gamma is the actual boundary velocity.
    beta = (interface_velocity_dot_n &
         -coeff%network_transport_fraction*vc_dot_n) &
         /coeff%network_diffusivity
  end function network_robin_beta

  pure elemental function free_robin_beta(coeff, vc_dot_n, &
       interface_velocity_dot_n) result(beta)
    type(actin_coefficients_t), intent(in) :: coeff
    real(dp), intent(in) :: vc_dot_n, interface_velocity_dot_n
    real(dp) :: beta

    beta = (interface_velocity_dot_n-vc_dot_n+coeff%membrane_rate) &
         /coeff%free_diffusivity
  end function free_robin_beta

  pure elemental function membrane_coupling_q(coeff) result(q)
    type(actin_coefficients_t), intent(in) :: coeff
    real(dp) :: q

    q = coeff%membrane_rate/coeff%network_diffusivity
  end function membrane_coupling_q

  pure elemental function network_robin_residual(coeff, theta_n, theta_c, &
       vc_dot_n, interface_velocity_dot_n, grad_theta_n_dot_n) result(residual)
    type(actin_coefficients_t), intent(in) :: coeff
    real(dp), intent(in) :: theta_n, theta_c, vc_dot_n
    real(dp), intent(in) :: interface_velocity_dot_n, grad_theta_n_dot_n
    real(dp) :: residual

    residual = grad_theta_n_dot_n &
         +network_robin_beta(coeff, vc_dot_n, interface_velocity_dot_n)*theta_n &
         -membrane_coupling_q(coeff)*theta_c
  end function network_robin_residual

  pure elemental function free_robin_residual(coeff, theta_c, vc_dot_n, &
       interface_velocity_dot_n, grad_theta_c_dot_n) result(residual)
    type(actin_coefficients_t), intent(in) :: coeff
    real(dp), intent(in) :: theta_c, vc_dot_n
    real(dp), intent(in) :: interface_velocity_dot_n, grad_theta_c_dot_n
    real(dp) :: residual

    residual = grad_theta_c_dot_n &
         +free_robin_beta(coeff, vc_dot_n, interface_velocity_dot_n)*theta_c
  end function free_robin_residual

  pure elemental function is_finite_coefficient(value) result(is_finite)
    real(dp), intent(in) :: value
    logical :: is_finite
    integer(int64) :: bits

    if (storage_size(value) == 64 .and. radix(value) == 2 .and. &
        digits(value) == 53 .and. minexponent(value) == -1021 .and. &
        maxexponent(value) == 1024) then
      bits = transfer(value, bits)
      is_finite = ibits(bits, 52, 11) /= int(z'7ff', int64)
    else
      is_finite = ieee_is_finite(value)
    endif
  end function is_finite_coefficient

end module actin_model_mod
