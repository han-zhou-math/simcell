! Implicit membrane--matrix adhesion on the moving cell boundary.
!
! The physical/current-arclength traction is
!
!   t_ad = -A_d V.
!
! newSpread integrates its marker argument over the fixed material coordinate
! ds, so the spread reference force is t_ad*J, J=|X_s|.  In contrast, the
! hydraulic slip law consumes the current traction t_ad and receives no J.
module fsi_adhesion_mod
  use parameters, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter, public :: ADHESION_OK = 0
  integer, parameter, public :: ADHESION_INVALID = 1
  public :: compute_adhesion_terms

contains

  subroutine compute_adhesion_terms(displacement_x,displacement_y,normal_x, &
       normal_y,jacobian,dt,adhesion,water_stress_mobility,reference_force_x, &
       reference_force_y,hydraulic_slip_x,hydraulic_slip_y,status)
    real(dp), intent(in) :: displacement_x(:),displacement_y(:)
    real(dp), intent(in) :: normal_x(:),normal_y(:),jacobian(:)
    real(dp), intent(in) :: dt,adhesion,water_stress_mobility
    real(dp), intent(out) :: reference_force_x(:),reference_force_y(:)
    real(dp), intent(out) :: hydraulic_slip_x(:),hydraulic_slip_y(:)
    integer, intent(out) :: status
    real(dp) :: velocity_x,velocity_y,normal_traction
    integer :: i,n

    status=ADHESION_INVALID
    reference_force_x=0.0_dp; reference_force_y=0.0_dp
    hydraulic_slip_x=0.0_dp; hydraulic_slip_y=0.0_dp
    n=size(displacement_x)
    if(size(displacement_y)/=n .or. size(normal_x)/=n .or. &
         size(normal_y)/=n .or. size(jacobian)/=n)return
    if(size(reference_force_x)/=n .or. size(reference_force_y)/=n)return
    if(size(hydraulic_slip_x)/=n .or. size(hydraulic_slip_y)/=n)return
    if(dt<=0.0_dp .or. adhesion<0.0_dp .or. water_stress_mobility<0.0_dp)return
    if(.not.ieee_is_finite(dt) .or. .not.ieee_is_finite(adhesion) .or. &
         .not.ieee_is_finite(water_stress_mobility))return
    if(any(jacobian<=0.0_dp))return
    if(.not.all(ieee_is_finite(displacement_x)) .or. &
         .not.all(ieee_is_finite(displacement_y)) .or. &
         .not.all(ieee_is_finite(normal_x)) .or. &
         .not.all(ieee_is_finite(normal_y)) .or. &
         .not.all(ieee_is_finite(jacobian)))return

    do i=1,n
      if(abs(normal_x(i)*normal_x(i)+normal_y(i)*normal_y(i)-1.0_dp)> &
           1.0e-8_dp)return
      velocity_x=displacement_x(i)/dt
      velocity_y=displacement_y(i)/dt
      ! Exactly one J converts the full current traction to a
      ! reference-coordinate force for newSpread.
      reference_force_x(i)=-adhesion*velocity_x*jacobian(i)
      reference_force_y(i)=-adhesion*velocity_y*jacobian(i)
      ! Water crosses only in the normal direction.  Its law consumes the
      ! normal component of current traction, hence deliberately no J here.
      normal_traction=-adhesion*(velocity_x*normal_x(i)+ &
           velocity_y*normal_y(i))
      hydraulic_slip_x(i)=water_stress_mobility*normal_traction*normal_x(i)
      hydraulic_slip_y(i)=water_stress_mobility*normal_traction*normal_y(i)
    enddo
    if(.not.all(ieee_is_finite(reference_force_x)) .or. &
         .not.all(ieee_is_finite(reference_force_y)) .or. &
         .not.all(ieee_is_finite(hydraulic_slip_x)) .or. &
         .not.all(ieee_is_finite(hydraulic_slip_y)))return
    status=ADHESION_OK
  end subroutine compute_adhesion_terms

end module fsi_adhesion_mod
