! SPDX-License-Identifier: BSD-3-Clause
!
! Source-compatible dualchem adapter for the neutral accepted-FSI snapshot.
!
! Dualchem's imported kernels still call the Stage-06 procedure names below.
! Stage 09 deliberately keeps that narrow source contract, but all stored data
! and interpolation now belong to fsi_transport_snapshot_mod.  This module owns
! no second velocity representation and performs no fluid or chemical solve.
module fsi_dualchem_velocity_bridge_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use parameters, only: dp,nx,ny,npts,zero,one
  use fsi_transport_snapshot_mod, only: fsi_transport_snapshot_t, &
      FSI_SNAPSHOT_OK,FSI_SNAPSHOT_NOT_READY,clear_fsi_transport_snapshot, &
      publish_fsi_transport_snapshot,copy_snapshot_face_velocity, &
      interpolate_snapshot_velocity
  implicit none
  private

  integer,parameter,public::VELOCITY_BRIDGE_OK=0
  integer,parameter,public::VELOCITY_BRIDGE_INVALID=1
  integer,parameter,public::VELOCITY_BRIDGE_NOT_READY=2

  type(fsi_transport_snapshot_t),save::installed_snapshot
  integer,save::installed_snapshot_id=0

  public::freeze_fsi_velocity
  public::copy_dualchem_face_velocity
  public::interpolate_fsi_velocity_to_marker
  public::clear_fsi_velocity_bridge
  public::install_dualchem_transport_snapshot
  public::get_installed_snapshot_id
  public::get_fsi_velocity_generation

contains

  subroutine clear_fsi_velocity_bridge()
    call clear_fsi_transport_snapshot(installed_snapshot)
    installed_snapshot_id=0
  end subroutine clear_fsi_velocity_bridge

  subroutine install_dualchem_transport_snapshot(snapshot, &
      expected_snapshot_id,status)
    ! Install one already validated accepted-FSI snapshot.  A failed or stale
    ! installation clears the adapter so no residual can reuse older velocity.
    type(fsi_transport_snapshot_t),intent(in)::snapshot
    integer,intent(in)::expected_snapshot_id
    integer,intent(out)::status
    real(dp),allocatable::vx_probe(:,:),vy_probe(:,:)
    integer::snapshot_status,allocation_status

    call clear_fsi_velocity_bridge()
    status=VELOCITY_BRIDGE_INVALID
    allocate(vx_probe(nx+1,ny),vy_probe(nx,ny+1),stat=allocation_status)
    if(allocation_status/=0)return
    call copy_snapshot_face_velocity(snapshot,vx_probe,vy_probe, &
        expected_snapshot_id,snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)return
    installed_snapshot=snapshot
    installed_snapshot_id=expected_snapshot_id
    status=VELOCITY_BRIDGE_OK
  end subroutine install_dualchem_transport_snapshot

  subroutine get_installed_snapshot_id(snapshot_id,status)
    integer,intent(out)::snapshot_id,status
    snapshot_id=0
    if(installed_snapshot_id<=0)then
      status=VELOCITY_BRIDGE_NOT_READY
      return
    end if
    snapshot_id=installed_snapshot_id
    status=VELOCITY_BRIDGE_OK
  end subroutine get_installed_snapshot_id

  subroutine get_fsi_velocity_generation(generation,status)
    ! Stage 6A moving-interface contract: the installed snapshot id is the
    ! velocity generation.  A trial that publishes a new accepted-FSI snapshot
    ! advances the generation; any affine-operator evaluation that observes a
    ! different generation is stale and must fail its inner solve.
    integer(int64),intent(out)::generation
    integer,intent(out)::status
    generation=0_int64
    if(installed_snapshot_id<=0)then
      status=VELOCITY_BRIDGE_NOT_READY
      return
    end if
    generation=int(installed_snapshot_id,int64)
    status=VELOCITY_BRIDGE_OK
  end subroutine get_fsi_velocity_generation

  subroutine freeze_fsi_velocity(u_mac,v_mac,status)
    ! Legacy/test entry point retained for the imported Stage-08 kernels and
    ! focused bridge tests.  Production Stage 09 installs a complete snapshot
    ! from the coupling harness instead.  Stationary placeholder geometry is
    ! sufficient here because the legacy API exposes velocity queries only.
    real(dp),intent(in)::u_mac(-1:,-1:),v_mac(-1:,-1:)
    integer,intent(out)::status
    type(fsi_transport_snapshot_t),allocatable::candidate
    real(dp)::x(npts),y(npts),normal(npts,2)
    integer::candidate_id,snapshot_status,allocation_status

    candidate_id=max(1,installed_snapshot_id+1)
    allocate(candidate,stat=allocation_status)
    if(allocation_status/=0)then
      call clear_fsi_velocity_bridge()
      status=VELOCITY_BRIDGE_INVALID
      return
    end if
    x=zero
    y=zero
    normal=zero
    normal(:,1)=one
    call publish_fsi_transport_snapshot(candidate,u_mac,v_mac,x,y,x,y,normal, &
        one,zero,candidate_id,snapshot_status)
    if(snapshot_status/=FSI_SNAPSHOT_OK)then
      call clear_fsi_velocity_bridge()
      status=VELOCITY_BRIDGE_INVALID
      return
    end if
    call install_dualchem_transport_snapshot(candidate,candidate_id,status)
  end subroutine freeze_fsi_velocity

  subroutine copy_dualchem_face_velocity(vx_face,vy_face,status)
    real(dp),intent(out)::vx_face(:,:),vy_face(:,:)
    integer,intent(out)::status
    integer::snapshot_status

    call copy_snapshot_face_velocity(installed_snapshot,vx_face,vy_face, &
        installed_snapshot_id,snapshot_status)
    call map_snapshot_status(snapshot_status,status)
  end subroutine copy_dualchem_face_velocity

  subroutine interpolate_fsi_velocity_to_marker(x,y,ux,uy,status)
    real(dp),intent(in)::x,y
    real(dp),intent(out)::ux,uy
    integer,intent(out)::status
    real(dp)::velocity(2)
    integer::snapshot_status

    call interpolate_snapshot_velocity(installed_snapshot,x,y,velocity, &
        installed_snapshot_id,snapshot_status)
    ux=velocity(1)
    uy=velocity(2)
    call map_snapshot_status(snapshot_status,status)
  end subroutine interpolate_fsi_velocity_to_marker

  subroutine map_snapshot_status(snapshot_status,bridge_status)
    integer,intent(in)::snapshot_status
    integer,intent(out)::bridge_status
    select case(snapshot_status)
    case(FSI_SNAPSHOT_OK)
      bridge_status=VELOCITY_BRIDGE_OK
    case(FSI_SNAPSHOT_NOT_READY)
      bridge_status=VELOCITY_BRIDGE_NOT_READY
    case default
      bridge_status=VELOCITY_BRIDGE_INVALID
    end select
  end subroutine map_snapshot_status

end module fsi_dualchem_velocity_bridge_mod
