! Stage 14 insertion of one frozen actin feedback object into the existing FSI
! right-hand side and membrane slip.
!
! newSpread remains the single immersed-boundary delta-kernel implementation.
! Its input is sigma*n*s_alpha; newSpread supplies material quadrature ds once.
! The interface slip is kw*sigma*n and therefore receives no Jacobian factor.
module fsi_actin_force_assembly_mod
  use parameters, only: dp,nx,ny,npts,xmin,xmax,ymin,ymax, &
       stage14_remove_uniform_actin_stress
  use IBmod, only: newSpread
  use scalar_operator_mod, only: is_finite_scalar
  use fsi_actin_feedback_mod, only: ACTIN_FEEDBACK_OK, &
       actin_fsi_feedback_t,validate_actin_fsi_feedback
  implicit none
  private

  integer,parameter,public::ACTIN_FORCE_OK=0
  integer,parameter,public::ACTIN_FORCE_INVALID=1
  public::add_actin_fsi_forcing

contains

  subroutine add_actin_fsi_forcing(feedback,marker_x,marker_y,normal_x, &
       normal_y,jacobian,water_permeability,rhs_u,rhs_v,slip_x,slip_y,status)
    type(actin_fsi_feedback_t),intent(in)::feedback
    real(dp),intent(in)::marker_x(:),marker_y(:),normal_x(:),normal_y(:)
    real(dp),intent(in)::jacobian(:),water_permeability
    real(dp),intent(inout)::rhs_u(-1:nx+1,-1:ny+1)
    real(dp),intent(inout)::rhs_v(-1:nx+1,-1:ny+1)
    real(dp),intent(out)::slip_x(:),slip_y(:)
    integer,intent(out)::status
    real(dp),allocatable::u_work(:,:),v_work(:,:),spread_u(:,:),spread_v(:,:)
    real(dp)::marker_force_x(npts),marker_force_y(npts)
    real(dp)::fluid_marker_stress(npts),uniform_stress
    real(dp)::slip_x_work(npts),slip_y_work(npts),normal_square,term,new_value
    integer::i,j,allocation_status
    logical::ok

    status=ACTIN_FORCE_INVALID
    slip_x=0.0_dp; slip_y=0.0_dp
    if(size(marker_x)/=npts .or. size(marker_y)/=npts)return
    if(size(normal_x)/=npts .or. size(normal_y)/=npts)return
    if(size(jacobian)/=npts)return
    if(size(slip_x)/=npts .or. size(slip_y)/=npts)return
    if(validate_actin_fsi_feedback(feedback)/=ACTIN_FEEDBACK_OK)return
    if(.not.all(is_finite_scalar(marker_x)))return
    if(.not.all(is_finite_scalar(marker_y)))return
    if(.not.all(is_finite_scalar(normal_x)))return
    if(.not.all(is_finite_scalar(normal_y)))return
    if(.not.all(is_finite_scalar(jacobian)))return
    if(.not.is_finite_scalar(water_permeability))return
    if(water_permeability<0.0_dp)return
    if(any(marker_x<=xmin) .or. any(marker_x>=xmax))return
    if(any(marker_y<=ymin) .or. any(marker_y>=ymax))return
    if(any(jacobian<=0.0_dp))return
    ! wrapLinSolve/TransRHS consume only the physical staggered unknowns.
    ! Ghost padding is solver workspace and can remain indeterminate, so it
    ! must not participate in validation or arithmetic here.
    if(.not.all(is_finite_scalar(rhs_u(0:nx-1,1:ny))))return
    if(.not.all(is_finite_scalar(rhs_v(1:nx,1:ny-1))))return

    allocate(u_work(-1:nx+1,-1:ny+1),v_work(-1:nx+1,-1:ny+1), &
         spread_u(-1:nx+1,-1:ny+1),spread_v(-1:nx+1,-1:ny+1), &
         stat=allocation_status)
    if(allocation_status/=0)return
    u_work=0.0_dp; v_work=0.0_dp
    u_work(0:nx-1,1:ny)=rhs_u(0:nx-1,1:ny)
    v_work(1:nx,1:ny-1)=rhs_v(1:nx,1:ny-1)
    marker_force_x=0.0_dp; marker_force_y=0.0_dp
    slip_x_work=0.0_dp; slip_y_work=0.0_dp
    fluid_marker_stress=feedback%marker_stress
    if(stage14_remove_uniform_actin_stress)then
      ! A constant isotropic stress changes only the pressure jump in
      ! incompressible Stokes flow.  Project it out of the regularized IB
      ! force to avoid an O(1/nu) parasitic velocity.  Keep the unprojected
      ! stress in the hydraulic slip below because water permeation responds
      ! to the absolute normal stress.
      uniform_stress=minval(fluid_marker_stress)
      fluid_marker_stress=fluid_marker_stress-uniform_stress
    endif

    do i=1,npts
      call checked_multiply(normal_x(i),normal_x(i),normal_square,ok)
      if(.not.ok)return
      call checked_multiply(normal_y(i),normal_y(i),term,ok)
      if(.not.ok)return
      call checked_add(normal_square,term,new_value,ok)
      if(.not.ok)return
      normal_square=new_value
      if(abs(normal_square-1.0_dp)>1.0e-8_dp)return

      call checked_multiply(fluid_marker_stress(i),normal_x(i),term,ok)
      if(.not.ok)return
      call checked_multiply(term,jacobian(i),marker_force_x(i),ok)
      if(.not.ok)return
      call checked_multiply(feedback%marker_stress(i),normal_x(i),term,ok)
      if(.not.ok)return
      call checked_multiply(water_permeability,term,slip_x_work(i),ok)
      if(.not.ok)return

      call checked_multiply(fluid_marker_stress(i),normal_y(i),term,ok)
      if(.not.ok)return
      call checked_multiply(term,jacobian(i),marker_force_y(i),ok)
      if(.not.ok)return
      call checked_multiply(feedback%marker_stress(i),normal_y(i),term,ok)
      if(.not.ok)return
      call checked_multiply(water_permeability,term,slip_y_work(i),ok)
      if(.not.ok)return
    enddo
    ! Keep every accumulation inside newSpread well away from overflow.  This
    ! bound is intentionally conservative and irrelevant at physical scales.
    if(maxval(abs(marker_force_x))>sqrt(huge(1.0_dp)))return
    if(maxval(abs(marker_force_y))>sqrt(huge(1.0_dp)))return
    call newSpread(marker_x,marker_y,marker_force_x,marker_force_y, &
         spread_u,spread_v)
    if(.not.all(is_finite_scalar(spread_u)))return
    if(.not.all(is_finite_scalar(spread_v)))return

    do j=1,ny
      do i=0,nx-1
        call checked_add(u_work(i,j),spread_u(i,j),new_value,ok)
        if(.not.ok)return
        u_work(i,j)=new_value
      enddo
    enddo
    do j=1,ny-1
      do i=1,nx
        call checked_add(v_work(i,j),spread_v(i,j),new_value,ok)
        if(.not.ok)return
        v_work(i,j)=new_value
      enddo
    enddo
    do j=1,ny
      do i=0,nx-1
        call checked_add(u_work(i,j),feedback%bulk_force_u(i+1,j), &
             new_value,ok)
        if(.not.ok)return
        u_work(i,j)=new_value
      enddo
    enddo
    do j=1,ny-1
      do i=1,nx
        call checked_add(v_work(i,j),feedback%bulk_force_v(i,j), &
             new_value,ok)
        if(.not.ok)return
        v_work(i,j)=new_value
      enddo
    enddo

    rhs_u(0:nx-1,1:ny)=u_work(0:nx-1,1:ny)
    rhs_v(1:nx,1:ny-1)=v_work(1:nx,1:ny-1)
    slip_x=slip_x_work; slip_y=slip_y_work
    status=ACTIN_FORCE_OK
  end subroutine add_actin_fsi_forcing

  subroutine checked_multiply(left,right,value,ok)
    real(dp),intent(in)::left,right
    real(dp),intent(out)::value
    logical,intent(out)::ok
    value=0.0_dp; ok=.false.
    if(.not.is_finite_scalar(left))return
    if(.not.is_finite_scalar(right))return
    if(left==0.0_dp .or. right==0.0_dp)then
      ok=.true.; return
    endif
    if(exponent(abs(left))+exponent(abs(right))>=maxexponent(value))return
    value=left*right
    ok=is_finite_scalar(value)
  end subroutine checked_multiply

  subroutine checked_add(left,right,value,ok)
    real(dp),intent(in)::left,right
    real(dp),intent(out)::value
    logical,intent(out)::ok
    real(dp)::sum
    value=0.0_dp; ok=.false.
    if(.not.is_finite_scalar(left))return
    if(.not.is_finite_scalar(right))return
    if(sign(1.0_dp,left)==sign(1.0_dp,right))then
      if(max(exponent(abs(left)),exponent(abs(right)))>= &
           maxexponent(value)-1)return
    endif
    sum=left+right
    if(.not.is_finite_scalar(sum))return
    value=sum; ok=.true.
  end subroutine checked_add

end module fsi_actin_force_assembly_mod
