! Topological side classification for a closed marker polygon.  The polygon is
! unwrapped once in the periodic x direction; each query is then shifted to the
! nearest periodic copy before applying an odd/even horizontal-ray test.
module interface_side_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use parameters, only: dp
  implicit none
  private

  integer,parameter,public::INTERFACE_SIDE_OK=0
  integer,parameter,public::INTERFACE_SIDE_INVALID=1
  public::prepare_periodic_polygon,classify_periodic_polygon_point
  public::classify_nearest_normal_point,stop_on_side_mismatch

contains

  subroutine prepare_periodic_polygon(marker_x,marker_y,period,polygon_x, &
       polygon_y,polygon_center,status)
    real(dp),intent(in)::marker_x(:),marker_y(:),period
    real(dp),intent(out)::polygon_x(:),polygon_y(:),polygon_center
    integer,intent(out)::status
    real(dp)::edge_dx
    integer::k,n

    status=INTERFACE_SIDE_INVALID
    polygon_x=0.0_dp
    polygon_y=0.0_dp
    polygon_center=0.0_dp
    n=size(marker_x)
    if(n<3 .or. size(marker_y)/=n)return
    if(size(polygon_x)/=n .or. size(polygon_y)/=n)return
    if(.not.ieee_is_finite(period).or.period<=0.0_dp)return
    if(.not.all(ieee_is_finite(marker_x)) .or. &
         .not.all(ieee_is_finite(marker_y)))return

    polygon_x(1)=marker_x(1)
    polygon_y=marker_y
    do k=2,n
      edge_dx=marker_x(k)-marker_x(k-1)
      edge_dx=edge_dx-period*anint(edge_dx/period)
      polygon_x(k)=polygon_x(k-1)+edge_dx
    enddo
    ! A physical cell is contractible in the periodic box.  After unwrapping,
    ! its closing edge must also be the nearest periodic connection.
    edge_dx=polygon_x(1)-polygon_x(n)
    if(abs(edge_dx-period*anint(edge_dx/period))>=0.5_dp*period)return
    polygon_center=sum(polygon_x)/real(n,dp)
    if(.not.ieee_is_finite(polygon_center))return
    status=INTERFACE_SIDE_OK
  end subroutine prepare_periodic_polygon

  subroutine classify_periodic_polygon_point(query_x,query_y,polygon_x, &
       polygon_y,period,polygon_center,is_inside,nearest_marker,distance_sq, &
       status)
    real(dp),intent(in)::query_x,query_y,polygon_x(:),polygon_y(:)
    real(dp),intent(in)::period,polygon_center
    logical,intent(out)::is_inside
    integer,intent(out)::nearest_marker
    real(dp),intent(out)::distance_sq
    integer,intent(out)::status
    real(dp)::shifted_x,trial_distance,crossing_x,denominator
    integer::k,previous,n

    status=INTERFACE_SIDE_INVALID
    is_inside=.false.
    nearest_marker=-1
    distance_sq=huge(1.0_dp)
    n=size(polygon_x)
    if(n<3 .or. size(polygon_y)/=n)return
    if(.not.ieee_is_finite(query_x).or. &
         .not.ieee_is_finite(query_y))return
    if(.not.ieee_is_finite(period).or.period<=0.0_dp)return
    ! The prepared polygon was validated once by prepare_periodic_polygon;
    ! avoid another O(nmarker) validation pass for every Eulerian query.
    if(.not.ieee_is_finite(polygon_center))return

    shifted_x=query_x+period*anint((polygon_center-query_x)/period)
    previous=n
    do k=1,n
      trial_distance=(shifted_x-polygon_x(k))**2+ &
           (query_y-polygon_y(k))**2
      if(trial_distance<=distance_sq)then
        distance_sq=trial_distance
        nearest_marker=k
      endif

      ! Half-open y intervals exclude horizontal edges and count a shared
      ! vertex exactly once, avoiding the usual ray/vertex ambiguity.
      if((polygon_y(k)>query_y).neqv.(polygon_y(previous)>query_y))then
        denominator=polygon_y(previous)-polygon_y(k)
        crossing_x=polygon_x(k)+(query_y-polygon_y(k))* &
             (polygon_x(previous)-polygon_x(k))/denominator
        if(shifted_x<crossing_x)is_inside=.not.is_inside
      endif
      previous=k
    enddo
    if(nearest_marker<1 .or. .not.ieee_is_finite(distance_sq))return
    status=INTERFACE_SIDE_OK
  end subroutine classify_periodic_polygon_point

  ! Reproduce the former side tag exactly: choose the nearest marker in raw
  ! Cartesian coordinates, then use the sign of the marker-normal projection.
  ! This is retained as a diagnostic reference while ray casting is active.
  subroutine classify_nearest_normal_point(query_x,query_y,marker_x,marker_y, &
       normal_x,normal_y,is_inside,nearest_marker,distance_sq,status)
    real(dp),intent(in)::query_x,query_y,marker_x(:),marker_y(:), &
         normal_x(:),normal_y(:)
    logical,intent(out)::is_inside
    integer,intent(out)::nearest_marker
    real(dp),intent(out)::distance_sq
    integer,intent(out)::status
    real(dp)::trial_distance,normal_projection
    integer::k,n

    status=INTERFACE_SIDE_INVALID
    is_inside=.false.
    nearest_marker=-1
    distance_sq=huge(1.0_dp)
    n=size(marker_x)
    if(n<1 .or. size(marker_y)/=n .or. size(normal_x)/=n .or. &
         size(normal_y)/=n)return
    if(.not.ieee_is_finite(query_x).or. &
         .not.ieee_is_finite(query_y))return
    if(.not.all(ieee_is_finite(marker_x)).or. &
         .not.all(ieee_is_finite(marker_y)).or. &
         .not.all(ieee_is_finite(normal_x)).or. &
         .not.all(ieee_is_finite(normal_y)))return

    do k=1,n
      trial_distance=(query_x-marker_x(k))**2+(query_y-marker_y(k))**2
      if(trial_distance<=distance_sq)then
        distance_sq=trial_distance
        nearest_marker=k
      endif
    enddo
    if(nearest_marker<1 .or. .not.ieee_is_finite(distance_sq))return
    normal_projection=(query_x-marker_x(nearest_marker))* &
         normal_x(nearest_marker)+(query_y-marker_y(nearest_marker))* &
         normal_y(nearest_marker)
    if(.not.ieee_is_finite(normal_projection))return
    is_inside=normal_projection<0.0_dp
    status=INTERFACE_SIDE_OK
  end subroutine classify_nearest_normal_point

  subroutine stop_on_side_mismatch(context,i,j,query_x,query_y,ray_inside, &
       nearest_inside,marker_x,marker_y)
    character(len=*),intent(in)::context
    integer,intent(in)::i,j
    real(dp),intent(in)::query_x,query_y,marker_x(:),marker_y(:)
    logical,intent(in)::ray_inside,nearest_inside
    integer::k,shape_unit,metadata_unit

    if(ray_inside.eqv.nearest_inside)return

    open(newunit=shape_unit,file='tagging_mismatch_shape.dat', &
         status='replace',action='write')
    write(shape_unit,'(a)')'# marker_x marker_y'
    do k=1,min(size(marker_x),size(marker_y))
      write(shape_unit,'(2(es24.16,1x))')marker_x(k),marker_y(k)
    enddo
    close(shape_unit)

    open(newunit=metadata_unit,file='tagging_mismatch.txt', &
         status='replace',action='write')
    write(metadata_unit,'(a)')'context='//trim(context)
    write(metadata_unit,'(a,i0)')'i=',i
    write(metadata_unit,'(a,i0)')'j=',j
    write(metadata_unit,'(a,es24.16)')'query_x=',query_x
    write(metadata_unit,'(a,es24.16)')'query_y=',query_y
    write(metadata_unit,'(a,l1)')'ray_inside=',ray_inside
    write(metadata_unit,'(a,l1)')'nearest_normal_inside=',nearest_inside
    close(metadata_unit)

    write(*,'(a)')'TAGGING_SIDE_MISMATCH'
    write(*,'(a)')'  context: '//trim(context)
    write(*,'(a,2(i0,1x))')'  grid indices: ',i,j
    write(*,'(a,2(es16.8,1x))')'  query position: ',query_x,query_y
    write(*,'(a,l1,a,l1)')'  ray inside: ',ray_inside, &
         '  nearest-normal inside: ',nearest_inside
    write(*,'(a)')'  membrane saved to tagging_mismatch_shape.dat'
    error stop 'Ray and nearest-normal side tags differ'
  end subroutine stop_on_side_mismatch

end module interface_side_mod
