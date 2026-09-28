! SPDX-License-Identifier: BSD-3-Clause
!
! Neutral, caller-owned transport state published after one FSI step accepts.
!
! The snapshot is a value object: publication copies the accepted MAC faces,
! marker geometry, and time metadata.  Chemical and actin consumers therefore
! read the same immutable time level and cannot alias or advance the FSI state.
module fsi_transport_snapshot_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use parameters, only: dp,nx,ny,npts,xmin,xmax,ymin,ymax,xlength,hg, &
      half,zero,one
  implicit none
  private

  integer,parameter,public::FSI_SNAPSHOT_OK=0
  integer,parameter,public::FSI_SNAPSHOT_INVALID=1
  integer,parameter,public::FSI_SNAPSHOT_NOT_READY=2
  integer,parameter,public::FSI_SNAPSHOT_STALE=3

  type,public::fsi_transport_snapshot_t
    private
    real(dp)::u_face(0:nx-1,1:ny)=zero
    real(dp)::v_face(1:nx,0:ny)=zero
    real(dp)::x_old(npts)=zero,y_old(npts)=zero
    real(dp)::x_new(npts)=zero,y_new(npts)=zero
    real(dp)::normal(npts,2)=zero
    real(dp)::interface_velocity(npts,2)=zero
    real(dp)::dt=zero,time=zero
    integer::snapshot_id=0
    logical::ready=.false.
  end type fsi_transport_snapshot_t

  public::publish_fsi_transport_snapshot
  public::copy_snapshot_face_velocity
  public::interpolate_snapshot_velocity
  public::copy_snapshot_marker_velocity
  public::copy_snapshot_interface_velocity
  public::copy_snapshot_marker_geometry
  public::clear_fsi_transport_snapshot

contains

  subroutine clear_fsi_transport_snapshot(snapshot)
    type(fsi_transport_snapshot_t),intent(inout)::snapshot
    snapshot%u_face=zero
    snapshot%v_face=zero
    snapshot%x_old=zero
    snapshot%y_old=zero
    snapshot%x_new=zero
    snapshot%y_new=zero
    snapshot%normal=zero
    snapshot%interface_velocity=zero
    snapshot%dt=zero
    snapshot%time=zero
    snapshot%snapshot_id=0
    snapshot%ready=.false.
  end subroutine clear_fsi_transport_snapshot

  subroutine publish_fsi_transport_snapshot(snapshot,u_mac,v_mac,x_old,y_old, &
      x_new,y_new,normal,dt,time,snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(inout)::snapshot
    real(dp),intent(in)::u_mac(-1:,-1:),v_mac(-1:,-1:)
    real(dp),intent(in)::x_old(:),y_old(:),x_new(:),y_new(:),normal(:,:)
    real(dp),intent(in)::dt,time
    integer,intent(in)::snapshot_id
    integer,intent(out)::status
    real(dp)::scale,tolerance,norm_sq(npts)

    call clear_fsi_transport_snapshot(snapshot)
    status=FSI_SNAPSHOT_INVALID
    if(size(u_mac,1)<nx+3 .or. size(u_mac,2)<ny+3)return
    if(size(v_mac,1)<nx+3 .or. size(v_mac,2)<ny+3)return
    if(size(x_old)/=npts .or. size(y_old)/=npts)return
    if(size(x_new)/=npts .or. size(y_new)/=npts)return
    if(size(normal,1)/=npts .or. size(normal,2)/=2)return
    if(snapshot_id<=0)return
    if(.not.is_finite_value(dt))return
    if(dt<=zero)return
    if(.not.is_finite_value(time))return
    if(.not.all(is_finite_value(u_mac(0:nx,0:ny+1))))return
    if(.not.all(is_finite_value(v_mac(0:nx,0:ny))))return
    if(.not.all(is_finite_value(x_old)))return
    if(.not.all(is_finite_value(y_old)))return
    if(.not.all(is_finite_value(x_new)))return
    if(.not.all(is_finite_value(y_new)))return
    if(.not.all(is_finite_value(normal)))return

    ! Markers must remain in the physical rectangle and in the interpolation
    ! support outside the half-cell wall band.  This is the same declared
    ! Stage-08 interface-away-from-wall limitation as the legacy bridge.
    if(any(x_old<xmin) .or. any(x_old>xmax))return
    if(any(x_new<xmin) .or. any(x_new>xmax))return
    if(any(y_old<ymin+half*hg) .or. any(y_old>ymax-half*hg))return
    if(any(y_new<ymin+half*hg) .or. any(y_new>ymax-half*hg))return

    norm_sq=normal(:,1)*normal(:,1)+normal(:,2)*normal(:,2)
    if(.not.all(is_finite_value(norm_sq)))return
    if(maxval(abs(norm_sq-one))>1024.0_dp*epsilon(one))return

    scale=max(one,maxval(abs(u_mac(0:nx,0:ny+1))), &
        maxval(abs(v_mac(0:nx,0:ny))))
    tolerance=64.0_dp*epsilon(one)*scale
    if(maxval(abs(u_mac(nx,0:ny+1)-u_mac(0,0:ny+1)))>tolerance)return
    if(maxval(abs(v_mac(0,0:ny)-v_mac(nx,0:ny)))>tolerance)return
    if(maxval(abs(v_mac(1:nx,0)))>tolerance)return
    if(maxval(abs(v_mac(1:nx,ny)))>tolerance)return
    if(maxval(abs(u_mac(0:nx-1,0)+u_mac(0:nx-1,1)))>tolerance)return
    if(maxval(abs(u_mac(0:nx-1,ny+1)+u_mac(0:nx-1,ny)))>tolerance)return

    snapshot%interface_velocity(:,1)=(x_new-x_old)/dt
    snapshot%interface_velocity(:,2)=(y_new-y_old)/dt
    if(.not.all(is_finite_value(snapshot%interface_velocity)))then
      call clear_fsi_transport_snapshot(snapshot)
      return
    end if

    snapshot%u_face=u_mac(0:nx-1,1:ny)
    snapshot%v_face=v_mac(1:nx,0:ny)
    snapshot%x_old=x_old
    snapshot%y_old=y_old
    snapshot%x_new=x_new
    snapshot%y_new=y_new
    snapshot%normal=normal
    snapshot%dt=dt
    snapshot%time=time
    snapshot%snapshot_id=snapshot_id
    snapshot%ready=.true.
    status=FSI_SNAPSHOT_OK
  end subroutine publish_fsi_transport_snapshot

  subroutine copy_snapshot_face_velocity(snapshot,vx_face,vy_face, &
      expected_snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    real(dp),intent(out)::vx_face(:,:),vy_face(:,:)
    integer,intent(in)::expected_snapshot_id
    integer,intent(out)::status

    vx_face=zero
    vy_face=zero
    status=FSI_SNAPSHOT_INVALID
    if(size(vx_face,1)/=nx+1 .or. size(vx_face,2)/=ny)return
    if(size(vy_face,1)/=nx .or. size(vy_face,2)/=ny+1)return
    call validate_snapshot_identity(snapshot,expected_snapshot_id,status)
    if(status/=FSI_SNAPSHOT_OK)return
    vx_face(1:nx,1:ny)=snapshot%u_face
    vx_face(nx+1,1:ny)=snapshot%u_face(0,1:ny)
    vy_face(1:nx,1:ny+1)=snapshot%v_face
  end subroutine copy_snapshot_face_velocity

  subroutine copy_snapshot_marker_velocity(snapshot,x,y,velocity,scale, &
      expected_snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    real(dp),intent(in)::x(:),y(:),scale
    real(dp),intent(out)::velocity(:,:)
    integer,intent(in)::expected_snapshot_id
    integer,intent(out)::status
    real(dp)::sample(2)
    integer::i,sample_status

    velocity=zero
    status=FSI_SNAPSHOT_INVALID
    if(size(x)/=npts .or. size(y)/=npts)return
    if(size(velocity,1)/=npts .or. size(velocity,2)/=2)return
    if(.not.all(is_finite_value(x)) .or. .not.all(is_finite_value(y)))return
    if(.not.is_finite_value(scale))return
    call validate_snapshot_identity(snapshot,expected_snapshot_id,status)
    if(status/=FSI_SNAPSHOT_OK)return
    do i=1,npts
      call interpolate_snapshot_velocity(snapshot,x(i),y(i),sample, &
          expected_snapshot_id,sample_status)
      if(sample_status/=FSI_SNAPSHOT_OK)then
        velocity=zero
        status=sample_status
        return
      end if
      velocity(i,:)=scale*sample
    end do
    if(.not.all(is_finite_value(velocity)))then
      velocity=zero
      status=FSI_SNAPSHOT_INVALID
      return
    end if
    status=FSI_SNAPSHOT_OK
  end subroutine copy_snapshot_marker_velocity

  subroutine copy_snapshot_interface_velocity(snapshot,velocity, &
      expected_snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    real(dp),intent(out)::velocity(:,:)
    integer,intent(in)::expected_snapshot_id
    integer,intent(out)::status

    velocity=zero
    status=FSI_SNAPSHOT_INVALID
    if(size(velocity,1)/=npts .or. size(velocity,2)/=2)return
    call validate_snapshot_identity(snapshot,expected_snapshot_id,status)
    if(status/=FSI_SNAPSHOT_OK)return
    velocity=snapshot%interface_velocity
  end subroutine copy_snapshot_interface_velocity

  subroutine copy_snapshot_marker_geometry(snapshot,x_new,y_new,normal, &
      expected_snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    real(dp),intent(out)::x_new(:),y_new(:),normal(:,:)
    integer,intent(in)::expected_snapshot_id
    integer,intent(out)::status

    ! Stage 10 needs the normals and marker positions from the same immutable
    ! value object as the face velocity.  Exposing a copy avoids a second,
    ! caller-maintained geometry channel and keeps snapshot components private.
    x_new=zero
    y_new=zero
    normal=zero
    status=FSI_SNAPSHOT_INVALID
    if(size(x_new)/=npts .or. size(y_new)/=npts)return
    if(size(normal,1)/=npts .or. size(normal,2)/=2)return
    call validate_snapshot_identity(snapshot,expected_snapshot_id,status)
    if(status/=FSI_SNAPSHOT_OK)return
    x_new=snapshot%x_new
    y_new=snapshot%y_new
    normal=snapshot%normal
  end subroutine copy_snapshot_marker_geometry

  subroutine interpolate_snapshot_velocity(snapshot,x,y,velocity, &
      expected_snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    real(dp),intent(in)::x,y
    real(dp),intent(out)::velocity(2)
    integer,intent(in)::expected_snapshot_id
    integer,intent(out)::status
    real(dp)::x_periodic,rx,ry,tx,ty
    integer::i0,i1,j0,j1,k0

    velocity=zero
    status=FSI_SNAPSHOT_INVALID
    call validate_snapshot_identity(snapshot,expected_snapshot_id,status)
    if(status/=FSI_SNAPSHOT_OK)return
    if(.not.is_finite_value(x))then
      status=FSI_SNAPSHOT_INVALID
      return
    end if
    if(.not.is_finite_value(y))then
      status=FSI_SNAPSHOT_INVALID
      return
    end if
    if(y<ymin+half*hg .or. y>ymax-half*hg)then
      status=FSI_SNAPSHOT_INVALID
      return
    end if

    x_periodic=xmin+modulo(x-xmin,xlength)
    rx=(x_periodic-xmin)/hg
    i0=floor(rx)
    if(i0==nx)i0=0
    tx=rx-real(i0,dp)
    i1=modulo(i0+1,nx)
    ry=(y-ymin)/hg-half
    if(ry>=real(ny-1,dp))then
      j0=ny-1
      ty=one
    else
      j0=floor(ry)+1
      ty=ry-real(j0-1,dp)
    end if
    j1=j0+1
    velocity(1)=(one-tx)*(one-ty)*snapshot%u_face(i0,j0)+ &
        tx*(one-ty)*snapshot%u_face(i1,j0)+ &
        (one-tx)*ty*snapshot%u_face(i0,j1)+ &
        tx*ty*snapshot%u_face(i1,j1)

    rx=(x_periodic-xmin)/hg-half
    k0=floor(rx)
    tx=rx-real(k0,dp)
    i0=modulo(k0,nx)+1
    i1=modulo(k0+1,nx)+1
    ry=(y-ymin)/hg
    j0=floor(ry)
    if(j0==ny)then
      j0=ny-1
      ty=one
    else
      ty=ry-real(j0,dp)
    end if
    j1=j0+1
    velocity(2)=(one-tx)*(one-ty)*snapshot%v_face(i0,j0)+ &
        tx*(one-ty)*snapshot%v_face(i1,j0)+ &
        (one-tx)*ty*snapshot%v_face(i0,j1)+ &
        tx*ty*snapshot%v_face(i1,j1)
    if(.not.all(is_finite_value(velocity)))then
      velocity=zero
      status=FSI_SNAPSHOT_INVALID
      return
    end if
    status=FSI_SNAPSHOT_OK
  end subroutine interpolate_snapshot_velocity

  subroutine validate_snapshot_identity(snapshot,expected_snapshot_id,status)
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    integer,intent(in)::expected_snapshot_id
    integer,intent(out)::status
    status=FSI_SNAPSHOT_NOT_READY
    if(.not.snapshot%ready)return
    status=FSI_SNAPSHOT_STALE
    if(expected_snapshot_id/=snapshot%snapshot_id)return
    status=FSI_SNAPSHOT_OK
  end subroutine validate_snapshot_identity

  elemental logical function is_finite_value(value)
    real(dp),intent(in)::value
    integer(int64)::bits,exponent_mask
    if(storage_size(value)==64)then
      bits=transfer(value,bits)
      exponent_mask=int(z'7FF0000000000000',int64)
      is_finite_value=iand(bits,exponent_mask)/=exponent_mask
    else
      is_finite_value=ieee_is_finite(value)
    end if
  end function is_finite_value

end module fsi_transport_snapshot_mod
