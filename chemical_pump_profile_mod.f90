! SPDX-License-Identifier: BSD-3-Clause
!
! Polarized active-solute pump used in the 2-D PNAS model.  The material
! coordinate is periodic on [0,2*pi): s=0 is the leading edge and s=pi the
! trailing edge.  A positive common amplitude therefore imports solute at the
! front (negative coefficient) and exports it at the rear (positive
! coefficient).  Equal Gaussian widths and amplitudes give the balanced PNAS
! profile from Eq. (III.17) of the supporting information.
module chemical_pump_profile_mod
  use parameters, only: dp, cpi, dualchem_pump_width, &
       dualchem_rear_width,dualchem_rear_amplitude_ratio
  implicit none
  private

  real(dp), parameter, public :: PNAS_PUMP_WIDTH=0.21_dp*cpi

  public :: pnas_gaussian_pump

contains

  pure elemental function pnas_gaussian_pump(s,common_amplitude) result(value)
    real(dp), intent(in) :: s, common_amplitude
    real(dp) :: value
    real(dp) :: front_profile,rear_profile
    integer :: image_index

    front_profile=exp(-(s/dualchem_pump_width)**2)+ &
        exp(-((s-2.0_dp*cpi)/dualchem_pump_width)**2)
    rear_profile=exp(-((s-cpi)/dualchem_pump_width)**2)
    if(dualchem_rear_width>0.0_dp)then
      rear_profile=0.0_dp
      do image_index=-2,2
        rear_profile=rear_profile+exp(-((s-cpi-2.0_dp*cpi*image_index)/dualchem_rear_width)**2)
      enddo
      rear_profile=dualchem_rear_amplitude_ratio*rear_profile
    endif
    value=common_amplitude*(-front_profile+rear_profile)
  end function pnas_gaussian_pump

end module chemical_pump_profile_mod
