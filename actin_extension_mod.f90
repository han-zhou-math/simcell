! Stage 10 scalar extension adapter for one actin species.
!
! This module deliberately reuses the accepted Stage 06 semi-periodic C++ MG
! entry point.  It does not introduce a second Cartesian solver or C ABI.  Its
! only operator transformation is the constant-diffusion identity
!
!   -D Lap(w)+div(v w)+lambda w=f
!       <=> -Lap(w)+div((v/D)w)+(lambda/D)w=f/D.
!
! The same normalization is used for the local correction equation, including
! dt_effective=D*dt and source_jump/D.  The physical inputs remain immutable;
! trial arrays are published only after correction, MG, and trace evaluation.
module actin_extension_mod
  use parameters, only: dp
  use grid_types, only: LagrangianGrid,EulerianGrid
  use scalar_operator_mod, only: scalar_operator_t,scalar_rhs_t, &
       SCALAR_STATUS_OK,validate_scalar_operator,validate_scalar_rhs, &
       compute_scalar_lambda,normalize_scalar_system,is_finite_scalar
  use geometry_mod, only: getCorrection,apply_full_correction, &
       apply_fresh_cleared_correction,valatIBpt
  use advec_diff_solver_mod, only: solve_advec_diff
  implicit none
  private

  integer,parameter,public::ACTIN_EXTENSION_OK=0
  integer,parameter,public::ACTIN_EXTENSION_INVALID=1
  integer,parameter,public::ACTIN_EXTENSION_BACKEND_FAILED=2

  type,public::interface_traces_t
    real(dp),allocatable::value_average(:),normal_average(:)
    real(dp),allocatable::value_interior(:),value_exterior(:)
    real(dp),allocatable::normal_interior(:),normal_exterior(:)
  end type interface_traces_t

  public::evaluate_actin_extension

contains

  subroutine evaluate_actin_extension(lag_geometry,tagged_geometry, &
       normal_jump,physical_operator,physical_rhs,dt,u_old,Cold, &
       trial_field,Cnew,crc,traces,status,backend_iterations)
    type(LagrangianGrid),intent(in)::lag_geometry
    type(EulerianGrid),intent(in)::tagged_geometry
    real(dp),intent(in)::normal_jump(:)
    type(scalar_operator_t),intent(in)::physical_operator
    type(scalar_rhs_t),intent(in)::physical_rhs
    real(dp),intent(in)::dt
    real(dp),intent(in)::u_old(:,:),Cold(:,:)
    real(dp),allocatable,intent(out)::trial_field(:,:),Cnew(:,:),crc(:,:)
    type(interface_traces_t),intent(out)::traces
    integer,intent(out)::status,backend_iterations

    type(scalar_operator_t)::op
    type(scalar_rhs_t)::rhs
    type(EulerianGrid)::work_grid
    real(dp),allocatable::field_work(:,:),Cwork(:,:),source(:,:),mask(:,:)
    real(dp),allocatable::zero_jump(:)
    real(dp)::lambda_physical,lambda_normalized,dt_effective
    real(dp)::inside(6),outside(6),dn_inside,dn_outside
    integer::nx_local,ny_local,nmarker,i,scalar_status,allocation_status
    logical::backend_converged

    status=ACTIN_EXTENSION_INVALID
    backend_iterations=0
    nx_local=tagged_geometry%nx_grid
    ny_local=tagged_geometry%ny_grid
    nmarker=lag_geometry%npts
    if(nx_local<1 .or. ny_local<1 .or. nmarker<1)return

    call allocate_zero_outputs(nx_local,ny_local,nmarker,trial_field,Cnew, &
         crc,traces,allocation_status)
    if(allocation_status/=0)return
    if(size(normal_jump)/=nmarker)return
    if(size(u_old,1)/=nx_local .or. size(u_old,2)/=ny_local)return
    if(size(Cold,1)/=nmarker .or. size(Cold,2)/=6)return
    if(.not.all(is_finite_scalar(normal_jump)))return
    if(.not.all(is_finite_scalar(u_old)))return
    if(.not.all(is_finite_scalar(Cold)))return
    if(validate_scalar_operator(physical_operator,nx_local,ny_local,nmarker) &
         /=SCALAR_STATUS_OK)return
    if(validate_scalar_rhs(physical_rhs,nx_local,ny_local,nmarker) &
         /=SCALAR_STATUS_OK)return

    call compute_scalar_lambda(physical_operator%reaction,.true.,dt, &
         lambda_physical,scalar_status)
    if(scalar_status/=SCALAR_STATUS_OK)return
    call normalize_scalar_system(physical_operator,physical_rhs, &
         lambda_physical,op,rhs,lambda_normalized,scalar_status)
    if(scalar_status/=SCALAR_STATUS_OK)return

    ! D*dt is the effective time scale of the unit-diffusion equation.  Split
    ! the guard so huge/dt is evaluated only when it cannot overflow.
    if(dt>1.0_dp)then
      if(physical_operator%diffusion>huge(1.0_dp)/dt)return
    endif
    dt_effective=physical_operator%diffusion*dt
    if(.not.is_finite_scalar(dt_effective) .or. dt_effective<=0.0_dp)return

    allocate(field_work(nx_local,ny_local),Cwork(nmarker,6), &
         source(nx_local,ny_local),mask(nx_local,ny_local), &
         zero_jump(nmarker),stat=allocation_status)
    if(allocation_status/=0)return
    field_work=0.0_dp; Cwork=0.0_dp; source=0.0_dp; mask=0.0_dp
    zero_jump=0.0_dp
    work_grid=tagged_geometry

    ! [theta]=0 and [dn theta]=psi.  The old correction is evaluated at the
    ! new marker about X^n by the inherited old-center repair in geometry_mod.
    call getCorrection(lag_geometry,zero_jump,normal_jump,rhs%source_jump, &
         dt_effective,0.0_dp,Cwork,Cold,op%reaction,.false.,.false., &
         op%velocity_marker)
    call apply_full_correction(lag_geometry,work_grid,Cwork,Cold,1,.false., &
         dt_effective)

    source=rhs%volume+u_old/dt_effective
    call apply_fresh_cleared_correction(lag_geometry,work_grid,Cold,source, &
         1,.false.,dt_effective,.true.,correct_both_sides=.true.)
    mask=0.5_dp*(real(work_grid%idf(1:nx_local,1:ny_local),dp)+1.0_dp)
    call solve_advec_diff(work_grid%x_min,work_grid%y_min,work_grid%x_max, &
         work_grid%y_max,nx_local,ny_local,op%vx_face,op%vy_face, &
         lambda_normalized,source,op%flux_n,op%flux_s,mask,work_grid%crc, &
         field_work,backend_iterations,backend_converged)
    if(backend_iterations<0)return
    ! Stage 6A backend contract: the multigrid solve declares convergence;
    ! a silent walk past the iteration ceiling is a backend failure here too.
    if(.not.backend_converged)then
      status=ACTIN_EXTENSION_BACKEND_FAILED
      return
    endif
    if(.not.all(is_finite_scalar(field_work)))then
      status=ACTIN_EXTENSION_BACKEND_FAILED
      return
    endif

    work_grid%u=0.0_dp
    work_grid%u(1:nx_local,1:ny_local)=field_work
    do i=1,nmarker
      call valatIBpt(inside,i,Cwork,lag_geometry,work_grid,work_grid%u,1)
      call valatIBpt(outside,i,Cwork,lag_geometry,work_grid,work_grid%u,-1)
      if(.not.all(is_finite_scalar(inside)) .or. &
           .not.all(is_finite_scalar(outside)))return
      traces%value_average(i)=0.5_dp*(inside(1)+outside(1))
      dn_inside=dot_product(lag_geometry%normal(i,:),inside(2:3))
      dn_outside=dot_product(lag_geometry%normal(i,:),outside(2:3))
      traces%normal_average(i)=0.5_dp*(dn_inside+dn_outside)
      traces%value_interior(i)=traces%value_average(i)
      traces%value_exterior(i)=traces%value_average(i)
      traces%normal_interior(i)=traces%normal_average(i)+0.5_dp*normal_jump(i)
      traces%normal_exterior(i)=traces%normal_average(i)-0.5_dp*normal_jump(i)
    enddo

    trial_field=field_work
    Cnew=Cwork
    crc=work_grid%crc
    status=ACTIN_EXTENSION_OK
  end subroutine evaluate_actin_extension

  subroutine allocate_zero_outputs(nx_local,ny_local,nmarker,field,Ccoef, &
       correction_grid,traces,status)
    integer,intent(in)::nx_local,ny_local,nmarker
    real(dp),allocatable,intent(out)::field(:,:),Ccoef(:,:),correction_grid(:,:)
    type(interface_traces_t),intent(out)::traces
    integer,intent(out)::status
    allocate(field(nx_local,ny_local),Ccoef(nmarker,6), &
         correction_grid(nx_local,ny_local),traces%value_average(nmarker), &
         traces%normal_average(nmarker),traces%value_interior(nmarker), &
         traces%value_exterior(nmarker),traces%normal_interior(nmarker), &
         traces%normal_exterior(nmarker),stat=status)
    if(status/=0)return
    field=0.0_dp; Ccoef=0.0_dp; correction_grid=0.0_dp
    traces%value_average=0.0_dp; traces%normal_average=0.0_dp
    traces%value_interior=0.0_dp; traces%value_exterior=0.0_dp
    traces%normal_interior=0.0_dp; traces%normal_exterior=0.0_dp
  end subroutine allocate_zero_outputs

end module actin_extension_mod
