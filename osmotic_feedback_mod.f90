! SPDX-License-Identifier: BSD-3-Clause
!
! Stage-07 physical osmotic feedback seam.
!
! This module contains only the chemical part of the membrane marker velocity
!
!   X_t,chem = kw * (c_i-c_e) * n,
!
! where n points from the interior to the exterior and kw already contains the
! model's RT scale.  Keeping this algebra in a PETSc-independent routine makes
! its sign and disabling limit directly testable; AdvanceFSI calls this exact
! routine before entering its implicit outer iteration.
module osmotic_feedback_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use parameters, only: dp
  implicit none
  private

  integer, parameter, public :: OSMOTIC_FEEDBACK_OK=0
  integer, parameter, public :: OSMOTIC_FEEDBACK_INVALID=1

  public :: compute_osmotic_slip_velocity

contains

  subroutine compute_osmotic_slip_velocity(chemical_jump,normal_x,normal_y, &
      kw_value,slip_x,slip_y,status)
    real(dp), intent(in) :: chemical_jump(:),normal_x(:),normal_y(:),kw_value
    real(dp), intent(out) :: slip_x(:),slip_y(:)
    integer, intent(out) :: status
    integer :: count

    slip_x=0.0_dp
    slip_y=0.0_dp
    status=OSMOTIC_FEEDBACK_INVALID
    count=size(chemical_jump)
    if(count<=0)return
    if(size(normal_x)/=count .or. size(normal_y)/=count .or. &
       size(slip_x)/=count .or. size(slip_y)/=count)return
    if(.not.ieee_is_finite(kw_value))return
    if(.not.all(ieee_is_finite(chemical_jump)))return
    if(.not.all(ieee_is_finite(normal_x)) .or. &
       .not.all(ieee_is_finite(normal_y)))return

    slip_x=kw_value*chemical_jump*normal_x
    slip_y=kw_value*chemical_jump*normal_y
    if(.not.all(ieee_is_finite(slip_x)) .or. &
       .not.all(ieee_is_finite(slip_y)))then
      slip_x=0.0_dp
      slip_y=0.0_dp
      return
    end if
    status=OSMOTIC_FEEDBACK_OK
  end subroutine compute_osmotic_slip_velocity

end module osmotic_feedback_mod
