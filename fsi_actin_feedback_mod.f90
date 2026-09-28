! Stage 14 accepted F-actin -> implicit-FSI feedback value object.
!
! This module is deliberately independent of PETSc, the immersed-boundary
! residual, and the scalar C++ backend.  It performs only the frozen mapping
!
!   accepted theta_n cell/marker state
!     -> MAC-face Brinkman drag alpha_B
!     -> explicit bulk force -eta D_n grad(theta_n)
!     -> marker network pressure sigma=k_sigma theta_n.
!
! The caller owns the accepted time level.  A successful object is therefore
! immutable by convention for one complete AdvanceFSI transaction.
module fsi_actin_feedback_mod
  use parameters, only: dp,nx,ny,npts,hg
  use scalar_operator_mod, only: is_finite_scalar
  use actin_model_mod, only: actin_coefficients_t
  implicit none
  private

  integer,parameter,public::ACTIN_FEEDBACK_OK=0
  integer,parameter,public::ACTIN_FEEDBACK_INVALID=1

  type,public::actin_fsi_feedback_t
    logical::enabled=.false.
    integer::source_snapshot_id=-1
    real(dp),allocatable::drag_u(:,:),drag_v(:,:)
    real(dp),allocatable::bulk_force_u(:,:),bulk_force_v(:,:)
    real(dp),allocatable::marker_stress(:)
  end type actin_fsi_feedback_t

  public::build_actin_fsi_feedback,build_disabled_actin_fsi_feedback
  public::validate_actin_fsi_feedback

contains

  subroutine build_disabled_actin_fsi_feedback(source_snapshot_id,feedback,status)
    integer,intent(in)::source_snapshot_id
    type(actin_fsi_feedback_t),intent(out)::feedback
    integer,intent(out)::status
    logical::initialized

    call initialize_zero_feedback(feedback,initialized)
    status=ACTIN_FEEDBACK_INVALID
    if(.not.initialized)return
    if(source_snapshot_id<0)return
    feedback%source_snapshot_id=source_snapshot_id
    status=ACTIN_FEEDBACK_OK
  end subroutine build_disabled_actin_fsi_feedback

  subroutine build_actin_fsi_feedback(auxiliary_field,correction_grid, &
       cell_side,u_face_side,v_face_side,marker_trace,coeff,eta,k_sigma, &
       source_snapshot_id,expected_source_snapshot_id,feedback,status)
    real(dp),intent(in)::auxiliary_field(:,:),correction_grid(:,:)
    integer,intent(in)::cell_side(:,:),u_face_side(:,:),v_face_side(:,:)
    real(dp),intent(in)::marker_trace(:)
    type(actin_coefficients_t),intent(in)::coeff
    real(dp),intent(in)::eta,k_sigma
    integer,intent(in)::source_snapshot_id,expected_source_snapshot_id
    type(actin_fsi_feedback_t),intent(out)::feedback
    integer,intent(out)::status
    real(dp),allocatable::interior_extension(:,:),drag_u_work(:,:), &
         drag_v_work(:,:),bulk_u_work(:,:),bulk_v_work(:,:),stress_work(:)
    real(dp)::drag_factor,bulk_factor
    real(dp)::theta_face,difference,gradient,value
    integer::i,j,left_cell,right_cell,allocation_status
    logical::ok,initialized

    call initialize_zero_feedback(feedback,initialized)
    status=ACTIN_FEEDBACK_INVALID
    if(.not.initialized)return

    ! Validate every scalar and array before relational tests or arithmetic;
    ! this ordering is required for signaling-NaN safety under the project's
    ! default invalid/zero/overflow traps.
    if(any(shape(auxiliary_field)/=[nx,ny]))return
    if(any(shape(correction_grid)/=[nx,ny]))return
    if(any(shape(cell_side)/=[nx,ny]))return
    if(any(shape(u_face_side)/=[nx,ny]))return
    if(any(shape(v_face_side)/=[nx,ny-1]))return
    if(size(marker_trace)/=npts)return
    if(.not.all(is_finite_scalar(auxiliary_field)))return
    if(.not.all(is_finite_scalar(correction_grid)))return
    if(.not.all(is_finite_scalar(marker_trace)))return
    if(.not.valid_coefficients(coeff))return
    if(.not.is_finite_scalar(eta))return
    if(.not.is_finite_scalar(k_sigma))return
    if(eta<=0.0_dp .or. k_sigma<=0.0_dp)return
    if(any(abs(cell_side)/=1))return
    ! These arrays are starter's native jdu/jdv masks, not dualchem side
    ! tags.  Their audited contract is binary: one means that the MAC face is
    ! in the physical interior and zero means that it is outside.  Keeping
    ! that convention here avoids a hidden data-layout translation in fmain.
    if(any(u_face_side<0) .or. any(u_face_side>1))return
    if(any(v_face_side<0) .or. any(v_face_side>1))return
    if(source_snapshot_id<0)return
    if(source_snapshot_id/=expected_source_snapshot_id)return

    allocate(interior_extension(nx,ny),drag_u_work(nx,ny), &
         drag_v_work(nx,ny-1),bulk_u_work(nx,ny),bulk_v_work(nx,ny-1), &
         stress_work(npts),stat=allocation_status)
    if(allocation_status/=0)return
    drag_u_work=0.0_dp; drag_v_work=0.0_dp
    bulk_u_work=0.0_dp; bulk_v_work=0.0_dp; stress_work=0.0_dp

    call checked_multiply(eta,1.0_dp-coeff%network_transport_fraction, &
         drag_factor,ok)
    if(.not.ok)return
    call checked_multiply(eta,coeff%network_diffusivity,bulk_factor,ok)
    if(.not.ok)return

    interior_extension=auxiliary_field
    do j=1,ny
      do i=1,nx
        if(cell_side(i,j)==-1)then
          call checked_add(auxiliary_field(i,j),correction_grid(i,j), &
               interior_extension(i,j),ok)
          if(.not.ok)return
        endif
      enddo
    enddo
    ! Negativity is an audit diagnostic, not an admissibility condition for
    ! the linear actin equation.  In particular, a smooth auxiliary extension
    ! may be negative away from the physical side without invalidating the PDE
    ! solve or the accepted time level.
    if(any(auxiliary_field<0.0_dp))then
      write(*,'(a,1x,es24.16)') 'STAGE14_NEGATIVE_AUXILIARY_EXTENSION', &
           minval(auxiliary_field)
    endif
    if(any(marker_trace<0.0_dp))then
      write(*,'(a,1x,es24.16)') 'STAGE14_NEGATIVE_NETWORK_TRACE', &
           minval(marker_trace)
    endif
    ! This clamp affects only the derived dissipative feedback coefficient;
    ! it does not modify or reject the accepted actin PDE state.  The Brinkman
    ! drag is required to remain nonnegative so it cannot become an artificial
    ! energy source when an auxiliary extension undershoots.
    where(interior_extension<0.0_dp)interior_extension=0.0_dp

    ! u-face storage has nx periodic faces.  Entry 1 is the seam between
    ! cell nx (left) and cell 1 (right); entry i>1 lies between i-1 and i.
    do j=1,ny
      do i=1,nx
        if(u_face_side(i,j)==0)cycle
        left_cell=modulo(i-2,nx)+1
        right_cell=modulo(i-1,nx)+1
        call checked_add(interior_extension(left_cell,j), &
             interior_extension(right_cell,j),value,ok)
        if(.not.ok)return
        call checked_multiply(0.5_dp,value,theta_face,ok)
        if(.not.ok)return
        call checked_add(interior_extension(right_cell,j), &
             -interior_extension(left_cell,j),difference,ok)
        if(.not.ok)return
        call checked_multiply(difference,1.0_dp/hg,gradient,ok)
        if(.not.ok)return
        call checked_multiply(drag_factor,theta_face, &
             drag_u_work(i,j),ok)
        if(.not.ok)return
        call checked_multiply(-bulk_factor,gradient, &
             bulk_u_work(i,j),ok)
        if(.not.ok)return
      enddo
    enddo

    ! v-face storage contains only the ny-1 physical faces between adjacent
    ! cell rows.  The solid walls themselves have no velocity unknown and no
    ! Stage-14 drag entry.
    do j=1,ny-1
      do i=1,nx
        if(v_face_side(i,j)==0)cycle
        call checked_add(interior_extension(i,j), &
             interior_extension(i,j+1),value,ok)
        if(.not.ok)return
        call checked_multiply(0.5_dp,value,theta_face,ok)
        if(.not.ok)return
        call checked_add(interior_extension(i,j+1), &
             -interior_extension(i,j),difference,ok)
        if(.not.ok)return
        call checked_multiply(difference,1.0_dp/hg,gradient,ok)
        if(.not.ok)return
        call checked_multiply(drag_factor,theta_face, &
             drag_v_work(i,j),ok)
        if(.not.ok)return
        call checked_multiply(-bulk_factor,gradient, &
             bulk_v_work(i,j),ok)
        if(.not.ok)return
      enddo
    enddo

    do i=1,npts
      call checked_multiply(k_sigma,max(marker_trace(i),0.0_dp), &
           stress_work(i),ok)
      if(.not.ok)return
    enddo

    ! Atomic publication boundary.  Every failure above leaves the initially
    ! allocated public object exactly zero.
    feedback%drag_u=drag_u_work
    feedback%drag_v=drag_v_work
    feedback%bulk_force_u=bulk_u_work
    feedback%bulk_force_v=bulk_v_work
    feedback%marker_stress=stress_work
    feedback%source_snapshot_id=source_snapshot_id
    feedback%enabled=.true.
    status=ACTIN_FEEDBACK_OK
  end subroutine build_actin_fsi_feedback

  subroutine initialize_zero_feedback(feedback,ok)
    type(actin_fsi_feedback_t),intent(out)::feedback
    logical,intent(out)::ok
    integer::allocation_status
    ok=.false.
    feedback%enabled=.false.
    feedback%source_snapshot_id=-1
    allocate(feedback%drag_u(nx,ny),feedback%drag_v(nx,ny-1), &
         feedback%bulk_force_u(nx,ny),feedback%bulk_force_v(nx,ny-1), &
         feedback%marker_stress(npts),stat=allocation_status)
    if(allocation_status/=0)return
    feedback%drag_u=0.0_dp
    feedback%drag_v=0.0_dp
    feedback%bulk_force_u=0.0_dp
    feedback%bulk_force_v=0.0_dp
    feedback%marker_stress=0.0_dp
    ok=.true.
  end subroutine initialize_zero_feedback

  integer function validate_actin_fsi_feedback(feedback) result(status)
    type(actin_fsi_feedback_t),intent(in)::feedback
    status=ACTIN_FEEDBACK_INVALID
    if(feedback%source_snapshot_id<0)return
    if(.not.allocated(feedback%drag_u) .or. &
         .not.allocated(feedback%drag_v))return
    if(.not.allocated(feedback%bulk_force_u) .or. &
         .not.allocated(feedback%bulk_force_v))return
    if(.not.allocated(feedback%marker_stress))return
    if(any(shape(feedback%drag_u)/=[nx,ny]))return
    if(any(shape(feedback%drag_v)/=[nx,ny-1]))return
    if(any(shape(feedback%bulk_force_u)/=[nx,ny]))return
    if(any(shape(feedback%bulk_force_v)/=[nx,ny-1]))return
    if(size(feedback%marker_stress)/=npts)return
    if(.not.all(is_finite_scalar(feedback%drag_u)))return
    if(.not.all(is_finite_scalar(feedback%drag_v)))return
    if(.not.all(is_finite_scalar(feedback%bulk_force_u)))return
    if(.not.all(is_finite_scalar(feedback%bulk_force_v)))return
    if(.not.all(is_finite_scalar(feedback%marker_stress)))return
    if(feedback%enabled)then
      if(any(feedback%drag_u<0.0_dp))return
      if(any(feedback%drag_v<0.0_dp))return
      if(any(feedback%marker_stress<0.0_dp))return
    else
      if(any(feedback%drag_u/=0.0_dp))return
      if(any(feedback%drag_v/=0.0_dp))return
      if(any(feedback%bulk_force_u/=0.0_dp))return
      if(any(feedback%bulk_force_v/=0.0_dp))return
      if(any(feedback%marker_stress/=0.0_dp))return
    endif
    status=ACTIN_FEEDBACK_OK
  end function validate_actin_fsi_feedback

  logical function valid_coefficients(coeff)
    type(actin_coefficients_t),intent(in)::coeff
    valid_coefficients=.false.
    if(.not.is_finite_scalar(coeff%network_transport_fraction))return
    if(.not.is_finite_scalar(coeff%network_diffusivity))return
    if(coeff%network_transport_fraction<=0.0_dp)return
    if(coeff%network_transport_fraction>=1.0_dp)return
    if(coeff%network_diffusivity<=0.0_dp)return
    valid_coefficients=.true.
  end function valid_coefficients

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
    ! This exponent bound is intentionally conservative by one endpoint: it
    ! rejects a possibly representable boundary product instead of allowing a
    ! rounded huge/right guard to pass and trap on overflow.
    if(exponent(abs(left))+exponent(abs(right))>=maxexponent(value))return
    value=left*right
    ok=is_finite_scalar(value)
  end subroutine checked_multiply

  subroutine checked_add(left,right,value,ok)
    real(dp),intent(in)::left,right
    real(dp),intent(out)::value
    logical,intent(out)::ok
    value=0.0_dp; ok=.false.
    if(.not.is_finite_scalar(left))return
    if(.not.is_finite_scalar(right))return
    if(sign(1.0_dp,left)==sign(1.0_dp,right))then
      if(max(exponent(abs(left)),exponent(abs(right)))>= &
           maxexponent(value)-1)return
    endif
    value=left+right
    ok=is_finite_scalar(value)
  end subroutine checked_add

end module fsi_actin_feedback_mod
