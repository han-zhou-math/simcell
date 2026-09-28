! SPDX-License-Identifier: BSD-3-Clause
!
! Balanced actin (de)polymerization profile from PNAS SI Eq. (III.12).
! The material coordinate is periodic on [0,2*pi): s=0 is the leading edge
! and s=pi is the trailing edge.  A positive amplitude polymerizes at the
! front and depolymerizes at the rear.
module actin_pnas_profile_mod
  use parameters, only: dp,cpi
  implicit none
  private

  real(dp),parameter,public::PNAS_ACTIN_WIDTH=0.21_dp*cpi

  public::pnas_actin_gaussian_profile

contains

  pure elemental function pnas_actin_gaussian_profile(s,amplitude) result(value)
    real(dp),intent(in)::s,amplitude
    real(dp)::value
    real(dp)::front_profile,rear_profile

    front_profile=exp(-(s/PNAS_ACTIN_WIDTH)**2)+ &
         exp(-((s-2.0_dp*cpi)/PNAS_ACTIN_WIDTH)**2)
    rear_profile=exp(-((s-cpi)/PNAS_ACTIN_WIDTH)**2)
    value=amplitude*(front_profile-rear_profile)
  end function pnas_actin_gaussian_profile

end module actin_pnas_profile_mod
