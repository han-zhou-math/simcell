! Distributed external load for two-dimensional load--velocity studies.
!
! The signed input external_load_x is the net x-directed force per unit
! out-of-plane depth, nondimensionalized by P0*L0.  The force is distributed
! uniformly per unit current arclength, so its integral around the membrane is
! exactly external_load_x.  A negative value opposes migration in +x.
module fsi_external_load_mod
  use parameters, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer,parameter,public::EXTERNAL_LOAD_OK=0
  integer,parameter,public::EXTERNAL_LOAD_INVALID=1
  public::compute_external_load_terms

contains

  subroutine compute_external_load_terms(normal_x,normal_y,jacobian, &
       material_spacing,external_load_x,water_stress_mobility, &
       reference_force_x,reference_force_y,hydraulic_slip_x, &
       hydraulic_slip_y,status)
    real(dp),intent(in)::normal_x(:),normal_y(:),jacobian(:)
    real(dp),intent(in)::material_spacing,external_load_x
    real(dp),intent(in)::water_stress_mobility
    real(dp),intent(out)::reference_force_x(:),reference_force_y(:)
    real(dp),intent(out)::hydraulic_slip_x(:),hydraulic_slip_y(:)
    integer,intent(out)::status
    real(dp)::perimeter,traction_x,normal_traction
    integer::i,n

    status=EXTERNAL_LOAD_INVALID
    reference_force_x=0.0_dp; reference_force_y=0.0_dp
    hydraulic_slip_x=0.0_dp; hydraulic_slip_y=0.0_dp
    n=size(jacobian)
    if(size(normal_x)/=n .or. size(normal_y)/=n)return
    if(size(reference_force_x)/=n .or. size(reference_force_y)/=n)return
    if(size(hydraulic_slip_x)/=n .or. size(hydraulic_slip_y)/=n)return
    if(n<3 .or. material_spacing<=0.0_dp .or. &
         water_stress_mobility<0.0_dp)return
    if(.not.ieee_is_finite(material_spacing) .or. &
         .not.ieee_is_finite(external_load_x) .or. &
         .not.ieee_is_finite(water_stress_mobility))return
    if(.not.all(ieee_is_finite(normal_x)) .or. &
         .not.all(ieee_is_finite(normal_y)) .or. &
         .not.all(ieee_is_finite(jacobian)))return
    if(any(jacobian<=0.0_dp))return

    perimeter=material_spacing*sum(jacobian)
    if(.not.ieee_is_finite(perimeter) .or. perimeter<=0.0_dp)return
    traction_x=external_load_x/perimeter
    do i=1,n
      if(abs(normal_x(i)*normal_x(i)+normal_y(i)*normal_y(i)-1.0_dp)> &
           1.0e-8_dp)return
      reference_force_x(i)=traction_x*jacobian(i)
      normal_traction=traction_x*normal_x(i)
      hydraulic_slip_x(i)=water_stress_mobility*normal_traction*normal_x(i)
      hydraulic_slip_y(i)=water_stress_mobility*normal_traction*normal_y(i)
    enddo
    if(.not.all(ieee_is_finite(reference_force_x)) .or. &
         .not.all(ieee_is_finite(hydraulic_slip_x)) .or. &
         .not.all(ieee_is_finite(hydraulic_slip_y)))return
    status=EXTERNAL_LOAD_OK
  end subroutine compute_external_load_terms

end module fsi_external_load_mod
