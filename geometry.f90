!---------------------------------------------------------------------------  
! module for setting up geometry for chemical module
! MODULE: main program 
!> @author 
!> Lingxing Yao
!
!DESCRIPTION: 
!@brief
!> Subroutines in this module handle geometric representation and grid crossing
!> of cell membrane with underling Eulerian grids. Those are needed for computing
!> chemical field
!> 
!> The subroutine getIBlist is for creating the linked list of locations for each
!> immersed boundary points (IB markers)
!--------------------------------------------------------------------------- 
module geometry
!!  USE ISO_C_BINDING
  use, intrinsic :: iso_c_binding
  use parameters
  use myfft
  use interface_side_mod, only: INTERFACE_SIDE_OK, &
      classify_nearest_normal_point
  implicit none

  private

  integer, parameter :: nxb = nx
  integer, parameter :: nyb = ny


  integer, dimension(1:nx*ny) :: idx
  integer, dimension(1:(nx-2)*(ny-2)) :: cip
  integer :: ctl, ctr, cbl, cbr !corners @ cell centers
  integer, dimension(ny-2) :: clf, crt
  integer, dimension(nx-2) :: cbt, ctp
!##
  type, public :: cc_bp
    type (cc_bp), pointer :: next, prev
    integer :: cor(1:4,1:2) !order lb, lt, rt, rb 
    integer :: icor(1:4)!indicating fluid(0)/cell(1) state
    integer :: id
!##    integer :: lbx, ltx, rbx, rtx, lby, lty, rby, rty
    double precision :: x, y ! physical coordinates of IB pt
    double precision :: vx, vy ! velocity of IB pt
    double precision :: ux, uy ! velocity of fluid @ IB pt
    double precision :: nx, ny ! normal @ current IB
    double precision :: err, sqr
    double precision :: pos(1:4,1:2)
    double precision :: vn(1:2) 
    double precision :: dummy
    double precision :: s
    integer :: k
    double precision :: ct, cts
    integer :: ij(6,2)
  end type cc_bp
!##
  type, public :: ibpt
    type (cc_bp), pointer :: p
  end type ibpt
!##

  type, public :: cc_augvar
    type (cc_augvar), pointer :: next, prev 
    double precision :: cm, cp, xb, yb, xm, ym, xp, yp, delm, delp, s, sz, sl, sm, sp, nx, ny
    double precision :: ux, uy, vdn,udn,dxs, fwF, ocm,ocp, ncm,ncp, nm,np, sgm, gsgm
    double precision :: st_m(4), st_p(4)
    integer :: ilocm(4,2),ilocp(4,2), k
    integer :: im(2), ip(2), iml, ipl!index of - & + side of (xb,yb)
    integer :: i_ali_grd_m(4,2), i_off_grd_m(4,2)
    integer :: i_ali_grd_p(4,2), i_off_grd_p(4,2)
    double precision :: p_ali_grd_m(3,2), p_off_grd_m(3,2)
    double precision :: p_ali_grd_p(3,2), p_off_grd_p(3,2)
    double precision :: c_ali_grd_m(4), c_off_grd_m(4)
    double precision :: c_ali_grd_p(4), c_off_grd_p(4)
    double precision :: stm_is(2), stm_js(2), stp_is(2), stp_js(2) ! for moving grid
    integer :: m_is(2), m_js(2), p_is(2), p_js(2) ! for moving grid
    integer :: mi_st, mj_st, pi_st, pj_st
  end type cc_augvar

  type, public :: iapt
    type (cc_augvar), pointer :: p
  end type iapt
!

  type, public :: mgtree
    double precision, pointer :: p(:,:)
    type(iapt), dimension(nent-2) :: aug
  end type mgtree


  public :: cip, clf, crt, cbt, ctp, cbl, cbr, ctl, ctr
  public :: getIBlist, BdyQuadParametric, curvatureInt, Tagging, FindBdyPtForFP,&
      getAugVarslist, paraval, smoothcrossing, getGeometry, releasemem, mydot, &
      dist, brs, sub2ind, cc_VelatIB, extend
contains
!!
subroutine getIBlist(ibary,dt)
! generate the linked IB list 
  type(cc_bp), pointer :: cb_old, curr
  type(ibpt), dimension(:) :: ibary(nent-2)
!
  integer :: il, ima, jma, i, j, info
  double precision :: x, y, dt, vx, vy, rh
!
  rh = one/hg

  !!do il = 1, nib !loop of identify platelet pts
  do il = 1, nent-2 !loop of identify platelet pts
    i = 2*nfil+1+(il-1)*nring
!=======================================================================
    allocate(cb_old) ! this initializes the 1st IB node & points to it
    cb_old%next => cb_old
    cb_old%prev => cb_old
    curr        => cb_old
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! We need to get the location of IB pts here, 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
    y = ypt(i); x = xpt(i)
    call cc_VelatIB(x, y, vx, vy, unc, vnc, info) ! interpolate fluid vel to IB pts
    curr%ux= vx; curr%uy= vy
    curr%x = x; curr%y = y
    curr%vx = (xpt(i)-oxpt(i))/dt; curr%vy = (ypt(i)-oypt(i))/dt
    ! Initial geometry is built before any spreading routine initializes the
    ! mutable xamin/yamin aliases.  Use the immutable domain origin directly:
    !     i=floor((x-xmin)/h+xshift),  j=floor((y-ymin)/h+yshift).
    ! On the original [0,1]^2 domain an uninitialized zero alias hid this bug;
    ! on the shared [-1,1]^2 domain it displaced every marker index by nx/2.
    ima = int((x-xmin)*rh+xshift); jma=int((y-ymin)*rh+yshift)
    cb_old%cor(1,1)=ima  ; cb_old%cor(1,2)=jma
    cb_old%cor(2,1)=ima  ; cb_old%cor(2,2)=jma+1
    cb_old%cor(3,1)=ima+1; cb_old%cor(3,2)=jma+1
    cb_old%cor(4,1)=ima+1; cb_old%cor(4,2)=jma
!
    do i = 2*nfil+2+(il-1)*nring, 2*nfil+nring+(il-1)*nring
      allocate(curr%next)
      y = ypt(i); x = xpt(i)
      call cc_VelatIB(x, y, vx, vy, unc, vnc, info)
      curr%next%ux= vx; curr%next%uy= vy
      curr%next%x = x; curr%next%y = y
      curr%next%vx = (xpt(i)-oxpt(i))/dt; curr%next%vy = (ypt(i)-oypt(i))/dt
      ima = int((x-xmin)*rh+xshift); jma=int((y-ymin)*rh+yshift)
!
      curr%next%cor(1,1)=ima  ; curr%next%cor(1,2) = jma
      curr%next%cor(2,1)=ima  ; curr%next%cor(2,2) = jma+1
      curr%next%cor(3,1)=ima+1; curr%next%cor(3,2) = jma+1
      curr%next%cor(4,1)=ima+1; curr%next%cor(4,2) = jma  
!##      print '(1x,"IDX: (",i3,",",i3,")")',ima, jma
      curr%next%prev => curr
      curr%next%next => cb_old
      curr           => curr%next
    enddo
!    print '(1x,"DIFF pos??",2(e13.6,1x)," ISEL",1x,i3)', tmp,tp1,isel
    curr%next%prev => curr ! this makes the list a doubly linked loop

!??  allocate(ibary(il)%p)
    !!ibary(il)%p => cb_old%prev
    !!ibary(il)%p => cb_old%next
    !!ibary(il)%p => cb_old
    ibary(il)%p => curr%next

!PPT    print *,ibary(il)%p%cor(1,:),curr%cor(1,:)
  enddo ! end of loop for allocating platelets pt

  return
end subroutine getIBlist
!
subroutine cc_VelatIB(x, y, vx, vy, un, vn, info)
! interpolate velocity field in (un,vn) to the position (x,y)
! and save it in (vx,vy)
  double precision :: x, y, vx, vy
  double precision, dimension(-1:nx+1,-1:ny+1) :: un, vn
  integer :: info, isel
!
  double precision :: tp1, tp2, dx, dy, rh, vcor(4,2), coe(4), xl, yl
  integer :: i, j, ik, jk
!
  dx = hg; dy = hg; rh = 1./hg;
!
!!    vx = 0.5*(u(i+1,j,1) + u(i,j,1));
!!    vy = 0.5*(u(i,j+1,2) + u(i,j,2))
!!    vmat(i,j,1) =-vx; vmat(i,j,2) =-vy
!!
  ik = int((x-xmin)*rh+xshift);   jk=int((y-ymin)*rh+yshift)
  xl = xmin+(dble(ik)-xshift)*dx; yl= ymin+(dble(jk)-yshift)*dy
  !!pos(1,1) = x;       pos(1,2) = y
  !!pos(2,1) = x;       pos(2,2) = y + dy
  !!pos(3,1) = x + dx;  pos(3,2) = y + dy
  !!pos(4,1) = x + dx;  pos(4,2) = y
!
  coe(1) = (xl+dx-x )*(yl+dy-y )
  coe(2) = (xl+dx-x )*(y    -yl)
  coe(3) = (x    -xl)*(y    -yl)
  coe(4) = (x    -xl)*(yl+dy-y )
!
  !!if (isel .eq. 0) then !x component
    vcor(1,1) = half*(un(ik  ,jk  ) + un(ik  ,jk+1))! (ik  ,jk  )
    vcor(2,1) = half*(un(ik  ,jk+1) + un(ik  ,jk+2))! (ik  ,jk+1)
    vcor(3,1) = half*(un(ik+1,jk+1) + un(ik+1,jk+2))! (ik+1,jk+1)
    vcor(4,1) = half*(un(ik+1,jk  ) + un(ik+1,jk+1))! (ik+1,jk  )

    vx = dot_product(coe,vcor(:,1))
    vx =-vx
  !!else if (isel .eq. 1) then ! y component
    vcor(1,2) = half*(vn(ik  ,jk  ) + vn(ik+1,jk  ))! (ik  ,jk  )
    vcor(2,2) = half*(vn(ik  ,jk+1) + vn(ik+1,jk+1))! (ik  ,jk+1)
    vcor(3,2) = half*(vn(ik+1,jk+1) + vn(ik+2,jk+1))! (ik+1,jk+1)
    vcor(4,2) = half*(vn(ik+1,jk  ) + vn(ik+2,jk  ))! (ik+1,jk  )
    vy = dot_product(coe,vcor(:,2))
    vy =-vy
  !!else
  !!  print *, 'should not be here, check calling options!'
  !!  stop
  !!endif
!
  return
end subroutine cc_VelatIB
!
!-------------------------------------------------------------------------------
!
subroutine BdyQuadParametric(iblst,mk,nb,uin)
!---+-B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
! To setup bdry representation by quadratic polynomial, following
! Balaras
!---+-B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
  implicit none
  type(ibpt), dimension(:) :: iblst(nent-2)
  double precision, intent(in out) :: mk(7,nring,nent-2), nb(nring,2, nent-2)
  double precision, intent(in out) :: uin(2,nring,nent-2)
!
  integer :: isel, iflg
  type(cc_bp), pointer :: curr, cb_bc
  double precision :: mI(3,3), mId(3,3), s(nring+1), mA(3,3),   &
          rx(3), ry(3), pc(3), pA(3,3), xx(3), rcond, mC(3), mR(3),     &
          ferr(1), berr(1), AF(3,3), work(24), rh
!!          ferr, berr, AF(3,3), work(24), rh
  integer :: ipiv(3), info, iwork(6), ij(6,2), kl(6,2), ir(6,2), strlen
  integer :: ima, jma, is, js, im, jm, cor(4,2), icor(4), ic(2)
  character :: equed
!
  integer :: jent, ip, jp, i1, i2
  integer i, j, k, il
  double precision :: tmp, tp1, tp2, sm, s0, s1, s2, ck, tupi, x0, y0, x1, y1, x2,  &
          y2 , x, y, ptc(4,2), alp, tp3, tp4, r0, r1, r2, r3, r4
!
  tupi = 2.d0*cpi
  rh = one/hg
!
  mI=0.d0; 
  do i = 1, 3
    mI(i,i) = 1.d0
  enddo
  mId=mI;
  do il = 1, nib
    curr => iblst(il)%p;
    sm = 0.d0; tmp = 0.d0; ck = 0.d0
!    do i = 1, nring
!      tmp = sqrt((curr%prev%x-curr%x)**2.d0+(curr%prev%y-curr%y)**2.d0)
!      sm  = sm + tmp
!      curr=> curr%prev
!    enddo
    !!prt print *, 'sm=== ', sm
    s(1) = 0.d0
    curr => iblst(il)%p;
    do i = 1, nring
      ! note we have the parameter s runs in (0,2pi)
      tmp = sqrt((curr%prev%x-curr%x)**2.d0+(curr%prev%y-curr%y)**2.d0)!*tupi/sm
      s(i+1) = s(i) + tmp
      ck = ck + tmp
      curr => curr%prev
    enddo
    s = s*tupi/s(nring+1)

    curr => iblst(il)%p;
    mA(:,1) = 1.d0
    do k = 1, nring
!      print *, 's(', k, ')', s(k)
      if (k .eq. 1) then
        s0 = s(nring) ; s1 = s(1) + tupi; s2 = s(2) + tupi;
        !!s0 = s(nring)-tupi ; s1 = s(1); s2 = s(2)
      else if (k .eq. nring) then
!        s0 = s(k-1); s1 = s(k); s2 = 0.d0
        s0 = s(k-1); s1 = s(k); s2 = tupi
      else
        s0 = s(k-1); s1 = s(k); s2 = s(k+1)
      endif
      x0 = curr%next%x; x1 = curr%x; x2 = curr%prev%x;
      y0 = curr%next%y; y1 = curr%y; y2 = curr%prev%y;
      mk(1,k,il) = s(k)
      curr%s  = s(k); curr%k = k
      mA(1,2) = s0; mA(1,3) = s0*s0;
      mA(2,2) = s1; mA(2,3) = s1*s1;
      mA(3,2) = s2; mA(3,3) = s2*s2;
      rx(1) = x0; rx(2) = x1; rx(3) = x2
      pc = rx
      pA = mA
      CALL DGESVX('Equilibration','No transpose',3,1,pA,3,AF,3,ipiv,   &
             equed,mR,mC,pc,3,xx,3,RCOND,FERR,BERR,WORK,IWORK,INFO)
      if (info .ne. 0) then
        print '(1x,12(e13.6,1x))', rcond
        stop
      endif
!
      pc = xx
      mk(2:4,k,il) = pc(1:3)
      ry(1) = y0; ry(2) = y1; ry(3) = y2
      pc = ry
      pA = mA
      CALL DGESVX('Equilibration','No transpose',3,1,pA,3,AF,3,ipiv,   &
             equed,mR,mC,pc,3,xx,3,RCOND,FERR,BERR,WORK,IWORK,INFO)
!!      call dgesv(3, 1, pA, 3, ipiv, pc, 3, info)
      if (info .ne. 0) then
        print '(1x,12(e13.6,1x))', rcond
        stop
      endif
!
      pc = xx
      mk(5:7,k,il) = pc(1:3)
      if ( k .eq. 1) then
        tp1 =  2.d0*mk(4,k,il)*tupi+mk(3,k,il)
        tp2 =  2.d0*mk(7,k,il)*tupi+mk(6,k,il)
      else
        tp1 =  2.d0*mk(4,k,il)*mk(1,k,il)+mk(3,k,il)
        tp2 =  2.d0*mk(7,k,il)*mk(1,k,il)+mk(6,k,il)
      endif
      tmp = sqrt(tp1*tp1+tp2*tp2)
      nb(k,1,il) = -tp2/tmp; nb(k,2,il) = tp1/tmp;
      mtd(k,1,il) =-tp2; mtd(k,2,il) = tp1;
      tdv(k,1,il)= -tp1/tmp; tdv(k,2,il)=-tp2/tmp
!=======================================================================
      curr%nx = -tp2/tmp; curr%ny = tp1/tmp;
      !!uin(1,k,il)=curr%vx-curr%ux;  uin(2,k,il)=curr%vy-curr%uy
      uin(1,k,il)=(-curr%vx+curr%ux)*curr%nx+(-curr%vy+curr%uy)*curr%ny
      !!uin(1,k,il)=-(curr%vx+curr%ux)*curr%nx -(curr%vy+curr%uy)*curr%ny
!=======================================================================
!      curr%vn(2)=nb(k,2,il)
      curr => curr%prev
    enddo
  enddo ! end of loop on different platelets for coef mk of polynomials
!========================================================================
! End of bdry quadratic polynomial representation
!========================================================================
  return
end subroutine BdyQuadParametric
!
!========================================================================
subroutine Tagging(ibary, id, idf, oid, kdf, mk, nb, isel)
! Tagging grid points in/out moving cell. 
! on input, oid & kdf should store id & idf from previous step;
! on input, mk & nb store bdry interpolation and normal respectively
! on input, isel = -1 for initial setup, when oid & kdf do not contain info
! on output, idf and id store current tagging info
!-----B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
  integer, intent(in out) :: id(-1:nxb+1,-1:nyb+1),  idf(-1:nxb+1,-1:nyb+1)
  integer, intent(in) :: oid(-1:nxb+1,-1:nyb+1), kdf(-1:nxb+1,-1:nyb+1)
  integer, intent(in) :: isel
  double precision, intent(in) :: mk(7,nring,nent-2), nb(nring,2,1:nent-2)
  type(ibpt), dimension(:) :: ibary(nent-2)
!
  logical :: is_intra, all_in(4)
  integer :: i,j,js,il,il_save,side_status,nearest_marker
  double precision :: x,y,dx,dy,tpv(nib)
  double precision :: nearest_distance_sq,old_marker_x(nring,nib), &
      old_marker_y(nring,nib)
  double precision :: ushift, vshift
  type(cc_bp),pointer::curr
!
!  print*, 'are we here in Tagging: mk1'
  dx = hg; dy = hg
  do il=1,nib
    ! The former implementation traversed the linked marker list through
    ! %prev.  Preserve that exact coordinate ordering when reproducing it;
    ! nb(:, :, il) was constructed in the same traversal order.
    curr=>ibary(il)%p
    do i=1,nring
      old_marker_x(i,il)=curr%x
      old_marker_y(i,il)=curr%y
      curr=>curr%prev
    enddo
  enddo
!!  if (time < dlt) then ! ANOTHER LAYER
  id = 0; idf = 0; iid = 0; IC_Tag = 0
    
  do j = 1, nyb
  do i = 1, nxb ! loop to identify intra/extra celullar pts over each cell center
    x = xmin+dx*(dble(i) - xshift); y = ymin+dy*(dble(j)-yshift)

    do il = 1, nib
      call classify_nearest_normal_point(x,y,old_marker_x(:,il), &
           old_marker_y(:,il),nb(:,1,il),nb(:,2,il),all_in(il), &
           nearest_marker,tpv(il),side_status)
      if(side_status/=INTERFACE_SIDE_OK)then
        error stop 'Nearest-normal side classification failed'
      endif
    enddo ! end of sweep on all IB objects and store in/out info in variable "all_in"

    il_save = minloc(tpv,1)
!
    is_intra = .false. 
    do il = 1, nib
!      if (al(il)) ijump = ijump + 1
      is_intra = is_intra .or. all_in(il)
      if (all_in(il)) IC_Tag(i,j) = il ! IMPORTANT! WE ASSUME GRID CENTERS ARE EXTRACELLULAR ON WITHIN 1 CELL MEMBRANE!
    enddo
!    if (ijump .ne. 0) then
    if (is_intra) then 
      ! if (i,j) is enclosed by ANY cell membrane idf(i,j) = 1, otherwise idf(i,j)=-1
!      id(i,j) = ijump
!!      id(i,j) = 1; idf(i,j) = 1
      idf(i,j) = 1
    else
!!      id(i,j) = 0; idf(i,j) =-1
      idf(i,j) =-1
      IC_Tag(i,j) =-il_save ! IMPORTANT! WE ASSUME CELL MEMBRANE DO NOT INTERSECT EACH OTHER!!
    endif
  enddo 
  enddo! end of loop for identify intra/extra cellular computational cell centers

!notneeded?  ! check carefully all the irregular points next internal interface
!notneeded?  do j = 2, nyb-1
!notneeded?  do i = 2, nxb-1 !IMPORTANT! WE ASSUME MOVING CELL IS 2 GRIDS AWAY FROM DOMAIN EDGE!
!notneeded?    js = abs(idf(i,j) - idf(i+1,j)) + abs(idf(i,j) - idf(i-1,j))        &
!notneeded?        +abs(idf(i,j) - idf(i,j+1)) + abs(idf(i,j) - idf(i,j-1))
!notneeded?    x = xmin+dx*(dble(i) - xshift); y = ymin+dy*(dble(j)-yshift)
!notneeded?    if (js .ne. 0) then
!notneeded?      do il = 1, nib
!notneeded?      !!print *, 'check now?', js, idf(i,j), i,j
!notneeded?        !call FindBdyPtForFP(s0,s1,il,x,y,xt,yt,mk,k,0)
!notneeded?        call FindBdyPtForFP(s0,s1,abs(il_save),x,y,xt,yt,mk,k,0)
!notneeded?
!notneeded?        tp1 = paraval(s1, nring, mk, il, 0, 1, info)
!notneeded?        tp2 = paraval(s1, nring, mk, il, 1, 1, info)
!notneeded?        tp4 = sqrt(tp1*tp1+tp2*tp2)
!notneeded?        vx  =-tp2; vy = tp1
!notneeded?        tmp = ((x-xt)*vx+(y-yt)*vy)/tp4/sqrt((x-xt)*(x-xt)+(y-yt)*(y-yt))
!notneeded?        tp3 = ((x-xt)*tp1+(y-yt)*tp2)/tp4/sqrt((x-xt)*(x-xt)+(y-yt)*(y-yt))
!notneeded?        ip = 0
!notneeded?        if (tmp <zero) then
!notneeded?          ip = 1
!notneeded?        else
!notneeded?          ip =-1
!notneeded?        endif
!notneeded?        if (tmp < zero ) then
!notneeded?!          id(i,j) = 1 ! inside of platelet
!notneeded?          all_in(il) = .true.      
!notneeded?        else
!notneeded?!          id(i,j) = 0
!notneeded?          all_in(il) = .false.
!notneeded?        endif
!notneeded?      enddo
!notneeded?!
!notneeded?      is_intra = .false. 
!notneeded?      ijump = 0
!notneeded?      do il = 1, nib
!notneeded?!        if (all_in(il)) ijump = ijump + 1
!notneeded?        is_intra = is_intra .or. all_in(il)
!notneeded?      enddo
!notneeded?!      if (ijump .ne. 0) then
!notneeded?      if (is_intra) then
!notneeded?!        id(i,j) = ijump
!notneeded?!!        id(i,j) = 1; idf(i,j) = 1
!notneeded?        ip = 1
!notneeded?      else
!notneeded?!!        id(i,j) = 0; idf(i,j) =-1
!notneeded?        ip =-1
!notneeded?      endif
!notneeded?      if (ip.ne. idf(i,j) ) then
!notneeded?        !!deb print '(1x,"not the same!", 4(i3,1x),1024(e16.9,1x))', i, j, &
!notneeded?        !!deb ip, idf(i,j), tmp, tp3, dist(x,y,xt,yt), dist(x,y,xt,yt)/hg, x, y, xt, yt, vx, vy
!notneeded?
!notneeded?        idf(i,j) = ip
!notneeded?      endif
!notneeded?    endif
!notneeded?  enddo
!notneeded?  enddo
!=======================================================================
!##
  id = 0
  do j = 1, nyb
  do i = 1, nxb! picking only grid pts just by the surface: forcing pt
!!    if (id(i,j) .eq. 1) then
!!    if (id(i,j) > 0 ) then
    js = abs(idf(i,j) - idf(i+1,j)) + abs(idf(i,j) - idf(i-1,j))        &
        +abs(idf(i,j) - idf(i,j+1)) + abs(idf(i,j) - idf(i,j-1))
    if (js .eq. 0) then
      id(i,j) = 0
    else
      id(i,j) = idf(i,j)
    endif
  enddo
  enddo! picking only grid pts just in the surface: forcing pts
  id(1,:)   = 0 ! this is 
  id(nxb,:) = 0 ! necessary 
  id(:,1)   = 0 ! to get
  id(:,nyb) = 0 ! first layer correct

!  iid(1:nx,1:ny) = idf(1:nx,1:ny)
  iid(1:nx,1:ny) = (idf(1:nx,1:ny)+1)/2
!  do j=1, ny
!  do i=1, nx
!    if(idf(i,j) < 0) then
!      iid(i,j) = 0 !iid set indicator for exterior points to be zero, interior to be one
!    endif
!  enddo
!  enddo
!!  ijd = iid
!!  do j = 1, ny
!!  do i = 1, nx
!!    if (id(i,j) .eq. -1) then
!!      ijd(i,j) = 1
!!    endif
!!  enddo
!!  enddo

  idu = 0; idv = 0; jdu = 0; jdv = 0
  ushift = zero; vshift = half
  do j = 2, ny-1
  do i = 1, nx-1
    if (dble(idf(i,j)+idf(i+1,j)) > -0.5) then
      idu(i,j) = 1
      if (dble(idf(i,j) + idf(i+1,j)) >1.1 .and. dble(idf(i,j+1)+idf(i+1,j+1))>1.1 .and. dble(idf(i,j-1)+idf(i+1,j-1))>1.1) then
        jdu(i,j) = 1
      endif
      x = xmin+dx*(dble(i) - ushift); y = ymin+dy*(dble(j)-vshift)

      do il = 1, nib
        call classify_nearest_normal_point(x,y,old_marker_x(:,il), &
             old_marker_y(:,il),nb(:,1,il),nb(:,2,il),all_in(il), &
             nearest_marker,nearest_distance_sq,side_status)
        if(side_status/=INTERFACE_SIDE_OK)then
          error stop 'Nearest-normal u-face classification failed'
        endif
      enddo ! end of sweep on platelets
!
      is_intra = .false. 
      do il = 1, nib
        is_intra = is_intra .or. all_in(il)
      enddo
      if (is_intra) then
        jdu(i,j) = 1
      else
        jdu(i,j) = 0
      endif
    endif
  enddo
  enddo
  ushift = half; vshift = 0
  do j = 2, ny-1
  do i = 2, nx-1
    if (idf(i,j) + idf(i,j+1) >-0.1) then
      idv(i,j) = 1
      if (dble(idf(i,j) + idf(i,j+1)) >1.1 .and. dble(idf(i+1,j)+idf(i+1,j+1))>1.1 .and. dble(idf(i-1,j)+idf(i-1,j+1))>1.1) then
        jdv(i,j) = 1
      endif
      x = xmin+dx*(dble(i) - ushift); y = ymin+dy*(dble(j)-vshift)

      do il = 1, nib
        call classify_nearest_normal_point(x,y,old_marker_x(:,il), &
             old_marker_y(:,il),nb(:,1,il),nb(:,2,il),all_in(il), &
             nearest_marker,nearest_distance_sq,side_status)
        if(side_status/=INTERFACE_SIDE_OK)then
          error stop 'Nearest-normal v-face classification failed'
        endif
      enddo ! end of sweep on platelets
!
      is_intra = .false. 
      do il = 1, nib
        is_intra = is_intra .or. all_in(il)
      enddo
      if (is_intra) then
        jdv(i,j) = 1
      else
        jdv(i,j) = 0
      endif
    endif
  enddo
  enddo
!========================================================================
! End of identifying fluid/solid points
!========================================================================
!  print *, 'Are we done in tagging: mk2'

  return
end subroutine Tagging
!
double precision function paraval(s, n, mk, il, idir, icase, info)
! evaluate function values along IB points, by using interpolation stored
! in mk(1:7,1:nring,1:nib), 
! il: the index of IB objects
! idir: 0<-x; 1<-y
! icase: 0<- original interpolation; 1<- 1st derivative; 2<- 2nd derivative
  integer :: lda, n
  double precision :: s, mk(7,n, nent-2)
  integer :: il, idir, icase, info
!
  double precision :: ck, tupi, eps
  integer :: i, k, inc, ijump
! 
  tupi = 2.*cpi
  select case (idir)
  case (0)
    inc=2
  case (1)
    inc=5
  case default
    print *, 'No such option, check carefully'
    stop
  end select
!
  eps = 1.d-10
  info = -1
  ! locate k
  ck = mod(mod(s, tupi)+tupi, tupi) ! archlength we will used
!  ck = s0
!PRT  print *, 's0', s0, 'ck', ck
  k = 1
  ijump = 0
  if ( ck > eps ) then
    if (ck < mk(1,k,il)) then
      ck = ck + tupi
      ijump = 1
      goto 654
    endif
    do while (k < n .and. ijump .eq. 0)
      !!if ( ck > mk(1,k,il) ) then 
      !!if ( ck > mk(1,k,il) ) then 
      !!  k = k + 1
      !!else  
      !!  ijump = 1
      !!endif  
      if ( ck .ge. mk(1,k,il) .and. ck < mk(1,k+1,il)) then
        ijump = 1
      else
        k = k + 1
      endif
    enddo
  else
    !!ck = ck + tupi
    k = 1
    ijump = 1
  endif
  if (k .eq. n) then 
    if (ijump .eq. 0) then
      if ( ck .ge. mk(1,k,il) .and. ck < tupi) then
        ijump = 1
      else
        print *, "double precisionly?", ck, s, mk(1,1,il)
        stop
      endif
    else
      print *, 'this could not happen!'
      stop
    endif
  endif
  if (k .eq. 1 .and. ijump .eq. 1) then
    !!print *, 'small ck', ck
    ck = tupi + ck
    !!ck = mod(ck, tupi)+tupi ! archlength we will used
  endif
  654 continue
  select case (icase)
  case (0)
    paraval = mk(inc,k,il) + (mk(inc+1,k,il) + mk(inc+2,k,il)*ck)*ck
  case (1)
!!  if (icase .eq. 1) then
    !!paraval = mk(inc,k,il) + mk(inc+1,k,il)*ck + mk(inc+2,k,il)*ck*ck
    paraval =  mk(inc+1,k,il) + 2.d0*mk(inc+2,k,il)*ck
!!  else if (icase .eq. -1) then ! linear interpolation
!!  endif
  !!prt print *, 'ck = ', ck, ' k = ', k, 'inc = ', inc, 'idir = ', idir, 's= ', s
  case (2)
    paraval = 2.d0*mk(inc+2,k,il)
  case default
    paraval = mk(inc,k,il) + (mk(inc+1,k,il) + mk(inc+2,k,il)*ck)*ck
    print *, 'no such option, evaluate original interpolation'
  end select
  if (ijump .eq. 1) then
    info = 0
  else
    print '(1x,"why here?", 4(i3,1x), 1024(e16.8,1x))',k, info, n, ijump, s, ck, paraval
    stop
  endif
!
  return
end function paraval
!
subroutine FindBdyPtForFP(s0,s1,il,x,y,xm,ym,mkk,k,isel)
!# isel may not be necessary, Jul-7-17
! we can not use ibary since previous storage are released
! on input (x,y) is the FP, and on output (xm,ym) is the
! corresponding bdry point for (x,y).
! on output s1 store the material parameter for (xm,ym)
implicit none
integer :: il, k,isel
double precision, intent(out) :: s1
double precision :: s0, x, y, xm, ym
!!type(ibpt), dimension(:) :: ibary(nent-2)
type(cc_bp), pointer :: curr, cb_bc, cb_lp

double precision, dimension(:,:,:) :: mkk(7,nring,nent-2)
double precision :: tmp, sk, xt, yt, x1, y1, x0, y0, tp1, tp2, ck,  x2, y2, tupi
integer :: kk, i, j, kkk, ip, lp, ijump, jump, info, nnt
double precision :: s2, s3, sl, su, so, sc, tp3, tp4, eps, erro
double precision :: opd, odis, nvl, nwl, nvu, nwu
  eps = 1.d-9
  erro = (1.d-12)*hg;
  !!erro = (1.d-13);
  nnt = 60
  jump = 0

  tupi = two*cpi
  tmp = 1.d10;
  kk = 1
!
!
  do i = 1, nring
    !!mk(1,i,il)
    x1 = paraval(mkk(1,i,il), nring, mkk,il,0,0,info)
    y1 = paraval(mkk(1,i,il), nring, mkk,il,1,0,info)
    tp1= ((x-x1)*(x-x1)+(y-y1)*(y-y1))
    if (tp1 <= tmp) then
      kk = i
      tmp=tp1
    endif
  enddo
  if (sqrt(tmp) > 1.44d0*hg) then
    print *, 'may be too agressive, can we handle?', kk, nring
    print '(1x,1024(e16.8,1x))', tmp, tmp/hg, x1, y1, x, y
  endif
  so = mkk(1,kk,il); odis = tmp
  if ( kk .eq. nring) then
    sl = mkk(1,kk-1,il); su = tupi+mkk(1,1,il);
  else if (kk .eq. 1) then
    sl = mkk(1,nring,il);su = tupi+mkk(1,2,il);
  else
    sl = mkk(1,kk-1,il); su = mkk(1,kk+1,il)
  endif
  !!print *, x, y, kk
!!  if (kk .eq. nring+1) then
!!    kk = nring
!!  else if (kk .eq. 0) then
!!    kk = 1
!!  endif
  kkk = kk
  s0 = mkk(1,kk,il); s1 = s0;
  !!xt = paraval(mkk(1,kk,il), nring, mkk,il,0,0,info); 
  !!yt = paraval(mkk(1,kk,il), nring, mkk,il,1,0,info); 

  call crs(nring,nib,il,mkk,s0,x1,y1,info)
  sk = s0 ! initial arclength
  kkk = info
!
  ip = 0
  402       continue
  ip=ip + 1
  tp1=f(mkk(4,kk,il),mkk(3,kk,il),mkk(2,kk,il),mkk(7,kk,il),                &
                  mkk(6,kk,il),mkk(5,kk,il),x,y,s0)
  tp2=fd(mkk(4,kk,il),mkk(3,kk,il),mkk(2,kk,il),mkk(7,kk,il),               &
                  mkk(6,kk,il),mkk(5,kk,il),x,y,s0)
  s2=s0
  s3=tp1/tp2
  s1=s0-tp1/tp2
!P          if (jp .eq. 1) print *, 's1', s1
!P          if (s1 < 0 .or. s1 > tupi) then
!P            print  '(1x,"s1",e14.6, 1x, 2(i3,1x), 10(e14.6,1x))', s1, ip,jp,&
!P                             tp1, tp2, tp1/tp2, mkk(1,jp,il)
!P            print *, 'ERROR!'
!P!            stop
!P          endif
  ck = mod(mod(s1, tupi)+tupi, tupi) ! archlength we will used
!  ck = s1 
  k  = 1
  ijump = 0
  if ( ck > eps ) then
    do while (k .le. nring .and. ijump .eq. 0)
      if ( ck > mkk(1,k,il) ) then 
        k = k + 1
      else  
        ijump = 1
      endif  
    enddo
  else
    ck = ck + tupi
    k = 1
  endif
  if (k .eq. nring+1) then 
!           print *, 'k = 161', ck
    k = 1
  endif
  s1 = ck
!102       continue
  call crs(nring,nib,il,mkk,s1,x2,y2,info)
  tp1 = dist(x1,y1,x2,y2)

  tp3 = paraval(s1, nring, mkk, il, 0, 1, info)
  tp4 = paraval(s1, nring, mkk, il, 1, 1, info)
  tp2 = sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2))
  tp1 = abs(((x-x2)/tp2*tp3+(y-y2)/tp2*tp4))/sqrt(tp3*tp3+tp4*tp4)

  if (tp1 > erro .and. ip < nnt) then
!          if (abs(s1-s0) > erro .and. ip < 20) then
!            call crs(nring,nib,il,mkk,ck,xm,ym,info)
!    tp1 = dist(xm,ym,x,y)
!     if (tp1 > 2.*h) then
!       s1 = s0 + pi
!       ip = ip + 1
!       goto 102
!     endif
    s0 = s1; kk = k; x1 = x2; y1 = y2
    goto 402
  endif
  !call crs(nring,nib,il,mkk,s1,xm,ym,k) ! find bdry pt for (i,j)
  xm = paraval(s1,nring,mkk,il,0,0, info)
  ym = paraval(s1,nring,mkk,il,1,0, info)

  if (dist(x,y,xm,ym) > 1.5*hg) then
    jump =  jump + 1
    !!print *,'Newton break reinitialize!'
    tp3=f(mkk(4,k ,il),mkk(3,k ,il),mkk(2,k ,il),mkk(7,k ,il),            &
                mkk(6,k ,il),mkk(5,k ,il),x,y,s0)
    tp4=f(mkk(4,k ,il),mkk(3,k ,il),mkk(2,k ,il),mkk(7,k ,il),            &
                mkk(6,k ,il),mkk(5,k ,il),x,y,s1)
    if (jump < 5) then
      ip = 0
      call random_number(tmp)
      s0 = mod(mod(sk + tmp*4.*cpi/dble(nring),tupi)+tupi,tupi)
      if (so < 1d-6) then
        sc = so
        call zoomin(sl,su,sc,x,y, il,kk,mkk,info)
        s0 = sc
        xt = paraval(sc, nring, mkk,il,0,0,info)
        yt = paraval(sc, nring, mkk,il,1,0,info)
        nvl= paraval(sc, nring, mkk,il,0,1,info)
        nwl= paraval(sc, nring, mkk,il,1,1,info)
        opd= nvl*(x-xt)+nwl*(y-yt)
        !!if (dabs(opd) < erro) then
        if (abs(opd) < eps) then
          s1 = s0; xm = xt; ym = yt
          goto 442
        endif
        !!prt print '(1x,"try", 1024(e16.8,1x))', s0, s1, mod(s0,tupi), x,y,xm,ym,dist(x,y,xm,ym)/h, opd
      endif
      goto 402
    endif

    print *, 'can not find the point', jump, isel
    !!print '(1x,1024(e15.8,1x))', sk, s0, s1, ck, tp1, tp2, tp3, tp4
    print '(1x,1024(e15.8,1x))', sk, s1, ck, x,y,xm,ym, dist(x,y,xm,ym)/hg
    print '(1x,1024(i3,1x))', kk, kkk, k, ip
    do il = 1, nib
    do i = 1, nring
      x1 = paraval(mkk(1,i,il), nring, mkk,il,0,0,info)
      y1 = paraval(mkk(1,i,il), nring, mkk,il,1,0,info)
      print '(1x,1024(e16.8,1x))', x1, y1
    enddo
    enddo
    info = -1
    stop
    goto 981
  endif
  442 continue
!
  nvl= paraval(s1, nring, mkk,1,0,1,info)
  nwl= paraval(s1, nring, mkk,1,1,1,info)
  opd=(nvl*(x-xm)+nwl*(y-ym))/sqrt(nvl*nvl+nwl*nwl)/sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym))
  if (abs(opd)> 1.d-3) then
    !!print *, 'might be something wrong with bdry pt, check it!'
    !!print *, opd, erro, s1, ip
    sc = so; info = 0
    call zoomin(sl,su,sc,x,y, il,kk,mkk,info)
    xt = paraval(sc, nring, mkk,il,0,0,info)
    yt = paraval(sc, nring, mkk,il,1,0,info)
    nvl= paraval(sc, nring, mkk,1,0,1,info)
    nwl= paraval(sc, nring, mkk,1,1,1,info)
    tmp=(nvl*(x-xt)+nwl*(y-yt))/sqrt(nvl*nvl+nwl*nwl)/sqrt((x-xt)*(x-xt)+(y-yt)*(y-yt))
    if (abs(tmp) < abs(opd)) then
      s1 = sc; xm = xt; ym = yt;
      !!print *, "corrected!", tmp, opd, erro, s1, ip
    else
      !!print '(1x,"hard to correct!  @", 1024(e16.8,1x))', tmp,opd,s1,x,y, xm,ym, dist(x,y,xm,ym)/hg, dble(ip)
    endif
    !stop
  endif
!
981 continue
  return
end subroutine FindBdyPtForFP
!
subroutine zoomin(sl,su,so,x,y, il,kk,mkk,info)
  double precision :: sl, su, so, x, y
  double precision, dimension(:,:,:) :: mkk(7,nring,nent-2)
  integer :: info, kk,il
!
  double precision :: xl, yl, xu, yu, vxl, vyl, vxu, vyu, vxc, vyc
  double precision :: tmp, sc, xc, yc, tpl, tpu, tpc, eps
  integer :: nmx, lp
!
  eps = 1d-7; nmx = 40; lp = 0
  xl = paraval(sl, nring, mkk,il,0,0,info)
  yl = paraval(sl, nring, mkk,il,1,0,info)
  vxl= paraval(sl, nring, mkk,il,0,1,info)
  vyl= paraval(sl, nring, mkk,il,1,1,info)
  xu = paraval(su, nring, mkk,il,0,0,info)
  yu = paraval(su, nring, mkk,il,1,0,info)
  vxu= paraval(su, nring, mkk,il,0,1,info)
  vyu= paraval(su, nring, mkk,il,1,1,info)
  
  tpl= vxl*(x-xl)+vyl*(y-yl)
  tpu= vxu*(x-xu)+vyu*(y-yu)
  if (tpl*tpu < zero) then
    lp = 0
    do while ( lp < nmx)
      sc = 0.5*(sl + su)
      xc = paraval(sc, nring, mkk,il,0,0,info)
      yc = paraval(sc, nring, mkk,il,1,0,info)
      vxc= paraval(sc, nring, mkk,il,0,1,info)
      vyc= paraval(sc, nring, mkk,il,1,1,info)
      tpc= vxc*(x-xc)+vyc*(y-yc)
      if (tpl*tpc <zero) then
        su = sc; xu = xc; yu = yc; tpu = tpc
      else
        sl = sc; xl = xc; yl = yc; tpl = tpc
      endif
      if (abs(su-sl) < eps) then
        so = sc
        exit
      else
        lp = lp + 1
      endif
    enddo
  else
    print *, 'not valid situation', tpl, tpu
    info = -1
  endif

!
  return
end subroutine zoomin
!
!
!1234567890123456789012345678901234567890123456789012345678901234567890/
subroutine crs(n,ni,id,mk,s0,x0,y0,k)
! find the pt (x0,y0) corresponds to arc length s0
! n: <-nring; ni:<-nib; id:<- il
! k: on output, return the segment number in mk where s0 is located
  implicit none
  integer :: n, ni, id
  double precision :: mk(7,n,ni),s0,x0,y0
!
  double precision :: ck, tupi, eps
  integer :: k, ijump

  eps = 1.d-8
!!  tupi=8.d0*atan(1.d0)
  tupi=2.d0*cpi
  ! locate k

  ck = mod(mod(s0, tupi)+tupi, tupi) ! archlength we will used
!  ck = s0
!PRT  print *, 's0', s0, 'ck', ck
  k = 1
  ijump = 0
  if ( ck > eps ) then
    do while (k .le. n .and. ijump .eq. 0)
      if ( ck > mk(1,k,id) ) then 
        k = k + 1
      else  
        ijump = 1
      endif  
    enddo
  else
    ck = ck + tupi
    k = 1
  endif
  if (k .eq. nring+1) then 
!     print *, 'k = 161', ck
    k = 1
  endif
  x0 = mk(2,k,id) + ck*(mk(3,k,id) + mk(4,k,id)*ck)
  y0 = mk(5,k,id) + ck*(mk(6,k,id) + mk(7,k,id)*ck)
!PRT  print *, 'x0', x0, 'y0', y0, 'k ', k 
!PRT  print '(1x, 3(e14.6,1x))', (mk(ijump,k,id),ijump=2,4)
!
  return
end subroutine crs

! 
double precision function f(ax,bx,cx,ay,by,cy,xi,yj,s)
  implicit none
  double precision :: ax,bx,cx,ay,by,cy,xi,yj,s
  
!!old  f = 2.d0*(ax*ax+ay*ay)*s*s*s + 3.d0*(ax*bx+ay*by)*s*s + (2.d0*ax*cx+  &
!!old      2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj))*s + (bx*cx+by*cy-bx*xi &
!!old      -by*yj)
  f = ((2.d0*(ax*ax+ay*ay)*s + 3.d0*(ax*bx+ay*by))*s + &
    (2.d0*ax*cx+2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj)))*s &
    + (bx*cx+by*cy-bx*xi -by*yj)
  return 
end function f

!

double precision function fd(ax,bx,cx,ay,by,cy,xi,yj,s)
! derivative of function f above
  implicit none
  double precision :: ax,bx,cx,ay,by,cy,xi,yj,s
  
!!old  fd = 6.d0*(ax*ax+ay*ay)*s*s + 6.d0*(ax*bx+ay*by)*s + (2.d0*ax*cx+     &
!!old      2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj))
  fd = (6.d0*(ax*ax+ay*ay)*s + 6.d0*(ax*bx+ay*by))*s + (2.d0*ax*cx+     &
      2.d0*ay*cy+bx*bx+by*by-2.d0*(ax*xi+ay*yj))
  return 
end function fd

double precision function dist(x1,y1,x2,y2)
  implicit none
  double precision :: x1, y1, x2, y2

  dist = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
!
  return
end function dist
!
!========================================================================
! compute the total curvature along cell
!========================================================================
subroutine curvatureInt(mk,conv, convsq)
  double precision, dimension(:,:,:) :: mk(7,nring,nent-2)
  double precision :: conv, convsq
!
  double precision :: xp, yp, xpp, ypp, tupi, s0, cv, dst, tcv, t1
  integer :: i, j, info, il
!
  tupi = 2.d0*cpi
  tcv = zero; 
  conv = zero
  convsq = zero
  dst = tupi/dble(nring)

  do il = 1, nib
    do i = 1, nring-1
      s0  = 0.5d0*(mk(1,i,il)+mk(1,i+1,il)); 
      dst = (mk(1,i+1,il)-mk(1,i,il))
      xp  = paraval(s0, nring, mk, il, 0, 1, info)
      yp  = paraval(s0, nring, mk, il, 1, 1, info)
      xpp = paraval(s0, nring, mk, il, 0, 2, info)
      ypp = paraval(s0, nring, mk, il, 1, 2, info)
      t1  = (xp*xp+yp*yp)
      cv  = (xp*ypp - yp*xpp)/t1*dst
      tcv = tcv + cv
      conv= conv + abs(cv)
      convsq = convsq + cv*cv
    enddo
    i = nring
    s0 = 0.5d0*(mk(1,i,il)+tupi); ! note, we have s0< tupi here
    dst = (tupi-mk(1,i,il))
    xp  = paraval(s0, nring, mk, il, 0, 1, info)
    yp  = paraval(s0, nring, mk, il, 1, 1, info)
    xpp = paraval(s0, nring, mk, il, 0, 2, info)
    ypp = paraval(s0, nring, mk, il, 1, 2, info)
    t1  = (xp*xp+yp*yp)
    cv  = (xp*ypp - yp*xpp)/t1*dst
    tcv = tcv + cv
    conv= conv + abs(cv)
    convsq = convsq + cv*cv
  enddo
  
!!  print '(1x,"absolute CURVATURE: ", e12.6, "curvature: ", 10(e12.6,1x))', conv/tupi, tcv/tupi, mk(1,1,il), mk(1,nring,il)
!
  return
end subroutine curvatureInt
!========================================================================
!
!========================================================================
subroutine getAugVarslist(ibary, iaary, idf, llen, mk, isel, info) ! genFP
  implicit none
  type(ibpt), dimension(:) :: ibary(nent-2)
  type(iapt), dimension(:) :: iaary(nent-2)
  !!double precision, dimension(:,:,:) :: uin(2,nring,nent-2)
  double precision, intent (in)  :: mk(7,nring,nent-2)
  integer, intent (in)  :: idf(-1:nxb+1,-1:nyb+1)
  integer :: isel, isd, ilq(nent-2)
  integer, intent (out) :: llen(nent-2)
  type(cc_bp), pointer :: cb_lp, cb_bc, curr
  type(cc_augvar), pointer :: cp_lp, cp_ss, cp_xx
!
!=======================================================================
  integer :: i, j, ix, iy, lp, lq, il, ik, jk, ip, jp, ict, is, js, idm
  integer :: jw, jq
  !!integer :: cor(4,2), icor(4), ic, info, lbx, lby, icase, jcase, iact, ks, ibt
  integer :: ic, info, lbx, lby, icase, jcase, iact, ks, ibt
  double precision :: tmp, tp1, tp2, s, sz, x, y, dx, dy, dell, delr, xi, yi, xt, yt
  double precision :: ptc(4,2), low, upp, nx, ny, dxs, smid, stmp, xmid, ymid
  double precision :: slow, supp, xlow, xupp, ylow, yupp, sc, s0, s1, xc, yc, sx, sy, sl
  double precision :: xo,yo, xp, yp, dtht, linkl
!
  dx = hg; dy =hg; dtht = two*cpi/dble(nring)
!
    !!deb il = 1
    !!deb !!curr => ibary(il)%p%next%next
    !!deb !!curr => ibary(il)%p%next
    !!deb curr => ibary(il)%p
    !!deb do i = 1, nring
    !!deb   print '(1x,4(2(e14.6,1x)))', curr%x, curr%y, xpt(i), ypt(i), curr%x-xpt(i), curr%y-ypt(i), curr%s
    !!deb   curr => curr%prev
    !!deb enddo

!deb  print *, 'Start finding GridX...'
  do il = 1, nib
    curr => ibary(il)%p%prev
    !!curr => ibary(il)%p%next
    !!curr => ibary(il)%p
    lq = 0; iact = 0
    lbx = curr%next%cor(1,1); lby = curr%next%cor(1,2)
    do i = 1, nring
      is = curr%cor(1,1) - lbx; js = curr%cor(1,2) - lby
      icase =abs(is)+abs(js) 
      slow = curr%next%s; supp = curr%s
      xlow = curr%next%x; xupp = curr%x
      ylow = curr%next%y; yupp = curr%y
      linkl= sqrt((xlow-xupp)*(xlow-xupp)+(ylow-yupp)*(ylow-yupp))

      !!deb print '(1x,1024(e12.6,1x))', slow, supp, xlow, ylow, xupp, yupp, curr%s
      !!deb if (i.eq. nring) stop
      !!prt print '(1x,"Dataset ", 1024(e16.8,1x))', slow, supp, xlow, xupp, ylow, yupp
      if (slow > supp .or. supp < 1d-8) then
        !!print '(1x,"we are here",2(i3,1x), 1024(e16.8,1x))', i,icase, slow, supp, xlow, ylow, xupp, yupp
        supp = 2.0d0*cpi-1d-13
        !!print '(1x,"after that ",2(i3,1x), 1024(e16.8,1x))', i,icase, slow, supp, xlow, ylow, xupp, yupp
        ks = nring;
        !!prt print '(1x,"supp= ", 1024(e15.8,1x))', slow, supp, 2.*cpi, mk(1,1,il), mk(1,nring, il)
        !!print *, 'here', i, mod(i+nring-2,nring),slow, supp
      endif
      !!if (i .eq. 2.or. i.eq. 3) then
      !!  print *, i, mod(i+nring-2,nring),slow, supp
      !!endif
      !!dxs = sqrt((xupp-xlow)**2.+(yupp-ylow)**2.)/abs(supp-slow)
      dxs = sqrt((xupp-xlow)*(xupp-xlow)+(yupp-ylow)*(yupp-ylow))/dtht
      select case (icase)
      case (0)
        !no grid crossing, do nothing but moving to next IB pt
      case (1)
        !!iact = iact + 1
        jcase = 2*is+js
        select case (jcase)
        case ( 2) ! (is,js)=(-1,0)
          xi = xmin+(dble(curr%cor(1,1))-xshift)*hg; 
          yi = ymin+(dble(curr%cor(1,2))-yshift)*hg
          xt = xi
          yt = ymin+(dble(curr%cor(2,2))-yshift)*hg
          low= xlow; upp= xupp; xc = xt; yc = yt
          ik = curr%cor(1,1); jk = curr%cor(1,2)
          ip = 0            ; jp = 1
          ic = 0
        case (-1) ! (is,js)=(0,1)
          xi = xmin+(dble(curr%cor(2,1))-xshift)*hg; 
          yi = ymin+(dble(curr%cor(2,2))-yshift)*hg
          xt = xmin+(dble(curr%cor(3,1))-xshift)*hg; 
          yt = yi
          low= ylow; upp= yupp; xc = xt; yc = yt
          ik = curr%cor(2,1); jk = curr%cor(2,2)
          ip = 1            ; jp = 0
          ic = 1
        case (-2) ! (is,js)=(1,0)
          xi = xmin+(dble(curr%cor(3,1))-xshift)*hg; 
          yi = ymin+(dble(curr%cor(3,2))-yshift)*hg
          xt = xi
          yt = ymin+(dble(curr%cor(4,2))-yshift)*hg
          low= xlow; upp= xupp; xc = xt; yc = yt
          ik = curr%cor(3,1); jk = curr%cor(3,2)
          ip = 0            ; jp =-1
          ic = 0
        case ( 1) ! (is,js)=(0,-1)
          xi = xmin+(dble(curr%cor(4,1))-xshift)*hg; 
          yi = ymin+(dble(curr%cor(4,2))-yshift)*hg
          xt = xmin+(dble(curr%cor(1,1))-xshift)*hg; 
          yt = yi
          low= ylow; upp= yupp; xc = xt; yc = yt
          ik = curr%cor(4,1); jk = curr%cor(4,2)
          ip =-1            ; jp = 0
          ic = 1
        case default
          print *, 'Not valid option here! stop'
          print '(1x,1024(e18.6,1x))', jcase, is, js, lbx, lby
          stop
        end select
        info = i
        call intercept(xc, yc, nx, ny, sx, sz, slow, supp, low, upp, mk, il, ic, info)
        sl = dtht*(dble(i-1)+sqrt((xlow-xc)*(xlow-xc)+(ylow-yc)*(ylow-yc))/linkl)
        !!sl = dtht*(dble(i-1)+(sx-slow)/dtht)
        !!prt print '(1x, "print A ", 1024(e16.8,1x))',  ylow, yc, yupp, xlow, xc, xupp, slow, sx, supp
        if (info .eq. -1) then
          print *, 'wrong return status, check 0', il, ic
          print '(1x,1024(e16.9,1x))', xc, yc, sx, slow,  supp, low, upp, xi,yi, xt, yt
          stop 
        endif
        ibt = 0
        if ( idf(ik,jk)* idf(ik+ip,jk+jp) > 0) then
          ibt = 1
          print '(1x,"special case: no grid points ", 1024(i3,1x))', idf(ik,jk), idf(ik+ip,jk+jp),  idf(ik,jk)* idf(ik+ip,jk+jp) 
        endif
        if (ibt .eq. 0) then
          if ((xt-xc)*(xc-xi)+(yt-yc)*(yc-yi)>0) then
            iact = iact + 1
            if (iact .eq. 1) then
              allocate(cp_ss); cp_xx => cp_ss
              cp_ss%k = i; cp_ss%s = sx; cp_ss%sz = sz
              call savegridcross(idf,xc, yc, nx, ny, sx, sl, dxs, ik, jk, ip, jp, cp_ss, info)
              call cc_errmsg(info, ibary, cp_ss)
            else
              allocate(cp_ss%next)
              cp_ss%next%k = i; cp_ss%next%s = sx; cp_ss%next%sz = sz
              cp_ss%next%prev=>cp_ss
              cp_ss%next%next=>cp_xx
              call savegridcross(idf,xc, yc, nx, ny, sx, sl, dxs, ik, jk, ip, jp, cp_ss%next, info)
              call cc_errmsg(info, ibary, cp_ss)
              cp_ss=>cp_ss%next
            endif
          else 
            print *, 'Very special case here, 3 grid crossing for one grid'
            print *, 'We just keep it here to see if it double precisionly happens'
            print *, ic, il, info, "info"
            print '(1x,1024(e16.8,1x))', xc, yc, xt, yt, xi, yi, slow, supp, sx
            print '(1x,1024(i3,1x))', is, js, icase, jcase, lbx, lby, idf(lbx,lby), ic
            print '(1x,1024(i3,1x))', (curr%cor(j,:),j=1,4), curr%icor
            print '(1x,1024(i3,1x))', (curr%next%cor(j,:),j=1,4), curr%next%icor
            curr => ibary(1)%p
            tmp = paraval(sx, nring, mk,il,1-ic,0,info);
            tp1 = paraval(slow, nring, mk,il,1-ic,0,info);
            tp2 = paraval(supp, nring, mk,il,1-ic,0,info);
            print '(1x,"now ", 1024(e16.9,1x))', tmp, sx, tp1, slow, tp2, supp, xc, yc, low, upp
            print '(1x,1024(e16.8,1x))', mk(1,1,il), mk(1,2,il), mk(1,3,il), mk(1,nring,il)
            print *, ' ', ibary(il)%p%s, ibary(il)%p%prev%s, ibary(il)%p%next%s
            print '(1x,1024(e16.8,1x))', curr%s, curr%prev%s, curr%next%s, curr%x, curr%y
            print *, ' '
            curr=>ibary(il)%p
            do j = 1, nring
              print '(1x,1024(e16.8,1x))', curr%x, curr%y, curr%s, mk(1, j,il)
              curr=>curr%prev
            enddo
            print '(1x,1024(e16.8,1x))', xc, yc, xt, yt, xi, yi, slow, supp, sx
            !!stop
          endif
!        else
!          print *, 'do not save this crossing!'
        endif


!=======================================================================
! special case for ik,jk to have 3 neighbor on other sides should be here
!=======================================================================
      case (2)
        ! 4 cases (is,js) = ( 1,-1):  
        !         (is,js) = (-1, 1): 
        !         (is,js) = ( 1, 1): 
        !         (is,js) = (-1,-1): 
        jcase = 2*is+js
        select case (jcase)
        case (-3)
          ik = curr%cor(3,1); jk = curr%cor(3,2); 
        case (-1)
          ik = curr%cor(4,1); jk = curr%cor(4,2); 
        case ( 3)
          ik = curr%cor(1,1); jk = curr%cor(1,2); 
        case ( 1)
          ik = curr%cor(2,1); jk = curr%cor(2,2); 
        case ( 2)
          if (js .eq. 2 .and. is .eq. 0) then
            print *, ''
            print '(1x,"show IB pt: ", 1024(e14.5,1x))', curr%x, curr%y
            do ik = 1, nring
              curr => curr%prev
              print '(1x,1024(e14.5,1x))', curr%x, curr%y
            enddo

            smid = 0.5*(slow+supp)
            xmid = paraval(smid, nring, mk,1,0,0,info)
            ymid = paraval(smid, nring, mk,1,1,0,info)

            print *, 'two IB pts are more than h apart'
            print '(1x,"icase & jcase ", 1024(i5,1x))', icase, jcase, is, js
            print '(1x, 1024(i5,1x))', lbx, lby, curr%cor(1,:)
            print '(1x, 1024(e14.6,1x))', slow, supp, xlow, xupp, ylow, yupp, xmid, ymid


!123456789212345678931234567894123456789512345678961234567897123456789
            stop
            
          else
            print '(1x,1024(e14.5,1x))', curr%x, curr%y
            do ik = 1, nring
              curr => curr%prev
              print '(1x,1024(e14.5,1x))', curr%x, curr%y
            enddo
            stop
          endif
        case default
          print *, 'It should not be here! stop', i, nring, iact
          !!print '(1x,1024(e18.6,1x))', jcase, is, js, lbx, lby
          print '(1x,1024(i6,1x))', jcase, is, js, lbx, lby
          print '(1x,1024(i6,1x))', curr%cor(1,:), curr%next%cor(1,:), curr%prev%cor(1,:)

          print '(1x,1024(e14.5,1x))', curr%x, curr%y

          do ik = 1, nring
            curr => curr%prev
            print '(1x,1024(e14.5,1x))', curr%x, curr%y
          enddo
          stop
        end select
        xt = xmin+(dble(ik)-xshift)*hg; xo = xt
        yt = ymin+(dble(jk)-yshift)*hg; yo = yt
        info = i
        call intercept(xt, yc, nx, ny, sx, sz, slow, supp, xlow, xupp, mk, il, 0, info)
        sl = dtht*(dble(i-1)+sqrt((xlow-xt)*(xlow-xt)+(ylow-yc)*(ylow-yc))/linkl)
        !!sl = dtht*(dble(i-1)+(sx-slow)/dtht)
        !!print '(1x, "print 1 ", 1024(e14.8,1x))',  curr%next%x, curr%next%y,  curr%x, curr%y, &
        !! curr%prev%x, curr%prev%y, curr%next%s, curr%s, curr%prev%s, slow, supp, xt, yt, dble(jcase)
        !!prt print '(1x, "print 1 ", 1024(e16.8,1x))',  ylow, yc, yupp, xlow, xt, xupp, slow, sx, supp
        if (info .eq. -1) then
          print *, 'wrong return status, check 1'
          stop 
        endif
        jp = int(dsign(1.d0,yc-yt)); ip = 0
        ibt = 0
        if ( idf(ik,jk)* idf(ik+ip,jk+jp) > 0) then
          ibt = 1 
	  iskip = .true.
          call FindBdyPtForFP(s0,s1,il,xo,yo,xp,yp,mk,info,isel)
          print '(1x,"special case: crossing skipped1", 1024(i5,1x))', idf(ik-ip,jk-jp), &
            idf(ik,jk), idf(ik+ip,jk+jp), ik, jk, ip, jp, curr%cor(1,:),  &
            curr%next%cor(1,:), lbx, lby, icase, jcase
!          print '(1x,1024(i3,1x))', ((idf(ik+jq,jk+jw),jq=-1,1),jw=-1,1)
!          print '(1x,1024(e16.9,1x))', xo-hg*dble(ip),yo-hg*dble(jp), xo, yo, xt, yc, &
!             xo+hg*dble(ip), yo+hg*dble(jp), xp, yp
          
          !!noprt cb_bc => curr
          !!noprt do info = 1, nring+1
          !!noprt   print '(1x,1024(e16.9,1x))', cb_bc%x, cb_bc%y
          !!noprt   cb_bc => cb_bc%prev
          !!noprt enddo
!          print '(1x,1024(e16.9,1x))', curr%prev%prev%prev%x,curr%prev%prev%prev%y,&
! curr%prev%prev%x, curr%prev%prev%y, curr%prev%x, curr%prev%y, curr%x, curr%y, &
! curr%next%x, curr%next%y, curr%next%next%x,  curr%next%next%y , curr%next%next%next%x,&
!  curr%next%next%next%y 
!          print '(1x,1024(e16.9,1x))', sx, slow, supp, xlow, ylow, xupp, yupp, dist(xo,yo,xt,yc)/hg
          !!deb curr => ibary(1)%p
          !!deb do lq = 1, nring
          !!deb   print '(1x,1024(e16.8,1x))', curr%x, curr%y
          !!deb   curr => curr%prev
          !!deb enddo
          !!deb stop
        endif
        if (ibt .eq. 0) then
          iact = iact + 1
          if (iact .eq. 1) then
            allocate(cp_ss); cp_xx => cp_ss
            cp_ss%k = i; cp_ss%s = sx; cp_ss%sz = sz
            call savegridcross(idf,xt, yc, nx, ny, sx, sl, dxs, ik, jk, ip, jp, cp_ss, info)
            call cc_errmsg(info, ibary, cp_ss)
          else
            allocate(cp_ss%next)
            cp_ss%next%k = i; cp_ss%next%s = sx; cp_ss%next%sz = sz
            cp_ss%next%prev=>cp_ss
            cp_ss%next%next=>cp_xx
            call savegridcross(idf,xt, yc, nx, ny, sx, sl, dxs, ik, jk, ip, jp, cp_ss%next, info)
            call cc_errmsg(info, ibary, cp_ss)
            cp_ss=>cp_ss%next
          endif
        !!else
        !!  print *, 'we skip this grid crossing!'
        endif

        !!if (ibt .eq. 1) then
        !!  print *, ' still need this?'
        !!endif
        info = i
        call intercept(xc, yt, nx, ny, sy, sz, slow, supp, ylow, yupp, mk, il, 1, info)
        sl = dtht*(dble(i-1)+sqrt((xlow-xc)*(xlow-xc)+(ylow-yt)*(ylow-yt))/linkl)
        !!sl = dtht*(dble(i-1)+(sy-slow)/dtht)
        !!print '(1x, "print 2 ", 1024(e14.8,1x))',  curr%next%x, curr%next%y,  curr%x, curr%y, &
        !! curr%prev%x, curr%prev%y, curr%next%s, curr%s, curr%prev%s, slow, supp, xt, yt, dble(jcase)
        !!prt print '(1x, "print 2 ", 1024(e16.8,1x))',  ylow, yt, yupp, xlow, xc, xupp, slow, sy, supp
        if (info .eq. -1) then
          print *, 'wrong return status, check 2'
          stop 
        endif
        jp = 0; ip = int(dsign(1.d0,dble(xc-xt)));
        ibt = 0
        if ( idf(ik,jk)* idf(ik+ip,jk+jp) > 0) then
          ibt = 1
	  iskip = .true.
          call FindBdyPtForFP(s0,s1,il,xo,yo,xp,yp,mk,info,isel)
          print '(1x,"special case: skipped here", 1024(i3,1x))', idf(ik-ip,jk-jp), &
            idf(ik,jk), idf(ik+ip,jk+jp), ik, jk, ip, jp, curr%cor(1,:),  &
            curr%next%cor(1,:), lbx, lby, icase, jcase, ik,jk,idn(ik,jk)
!          print '(1x,1024(i3,1x))', ((idf(ik+jq,jk+jw),jq=-1,1),jw=-1,1)
!          print '(1x,1024(e16.9,1x))', xo-hg*dble(ip),yo-hg*dble(jp), xo, yo, xc, yt, &
!            xo+hg*dble(ip), yo+hg*dble(jp), xp, yp

          !!noprt cb_bc => curr
          !!noprt do info = 1, nring+1
          !!noprt print '(1x,1024(e16.9,1x))', cb_bc%x, cb_bc%y
          !!noprt cb_bc => cb_bc%prev
          !!noprt enddo
          !!print '(1x,1024(e16.9,1x))', curr%prev%prev%prev%x,curr%prev%prev%prev%y,  curr%prev%prev%x, curr%prev%prev%y, curr%prev%x, curr%prev%y, curr%x, curr%y, curr%next%x, curr%next%y, curr%next%next%x,  curr%next%next%y , curr%next%next%next%x, curr%next%next%next%y 
!          print '(1x,1024(e16.9,1x))', sy, slow, supp, xlow, ylow, xupp, yupp, dist(xo,yo,xc,yt)/hg
        endif
        if (ibt .eq. 0) then
          iact = iact + 1
          if (iact .eq. 1) then
            allocate(cp_ss); cp_xx => cp_ss
            cp_ss%k = i; cp_ss%s = sy; cp_ss%sz = sz
            call savegridcross(idf,xc, yt, nx, ny, sy, sl, dxs, ik, jk, ip, jp, cp_ss, info)
            call cc_errmsg(info, ibary, cp_ss)
          else
            allocate(cp_ss%next)
            cp_ss%next%k = i; cp_ss%next%s = sy; cp_ss%next%sz = sz
            cp_ss%next%prev=>cp_ss
            cp_ss%next%next=>cp_xx
            call savegridcross(idf,xc, yt, nx, ny, sy, sl, dxs, ik, jk, ip, jp, cp_ss%next, info)
            call cc_errmsg(info, ibary, cp_ss)
            cp_ss=>cp_ss%next
          endif
          !!old allocate(cp_ss%next)
          !!old cp_ss%next%prev=>cp_ss
          !!old cp_ss%next%next=>cp_xx
          !!old cp_ss%next%k = i; cp_ss%next%s = sy
          !!old call savegridcross(xc, yt, nx, ny, sy, ik, jk, ip, jp, cp_ss%next, info)
          !!old call cc_errmsg(info, ibary, cp_ss)
          !!old cp_ss => cp_ss%next
        !!else
        !!  print *, 'we skip this here!'
        endif
    687 continue

      case default
        print *, 'Wrong cases, there is something missing'
        print '(1x,1024(i5,1x))', curr%k, i
        print '(1x,1024(i5,1x))', is, js, lbx, lby, curr%cor(1,:), icase, jcase
        print '(1x,1024(i5,1x))', ((curr%cor(is,js),js=1,2),is=1,4)
        print '(1x,1024(i5,1x))', ((curr%next%cor(is,js),js=1,2),is=1,4)
        print '(1x,1024(i5,1x))', ((curr%next%next%cor(is,js),js=1,2),is=1,4)
        print *, ' '
        curr =>ibary(il)%p
        do is = 1, nring
          print '(1x, 12(i3,1x),2(e13.6,1x))', is, curr%k, lbx, lby, (curr%cor(j,:),j=1,4), curr%x, curr%y
          curr => curr%prev
        enddo
        stop
      end select 
      lbx = curr%cor(1,1); lby = curr%cor(1,2)
      curr =>curr%prev
      !!curr =>curr%next
    enddo
!!    print *, '# of grid crossing', iact
    cp_ss%next%prev => cp_ss

    !!iaary(il)%p => cp_ss%next
    iaary(il)%p => cp_ss

    llen(il) = iact
    !!deb print '(1x,"pre   ", 4(i3,1x), 1024(e16.8,1x))', iaary(il)%p%prev%im, iaary(il)%p%prev%ip,iaary(il)%p%prev%xb,iaary(il)%p%prev%yb,iaary(il)%p%prev%s
    !!deb print '(1x,"First ", 4(i3,1x), 1024(e16.8,1x))', iaary(il)%p%im, iaary(il)%p%ip,iaary(il)%p%xb,iaary(il)%p%yb,iaary(il)%p%s
    !!deb print '(1x,"next  ", 4(i3,1x), 1024(e16.8,1x))', iaary(il)%p%next%im, iaary(il)%p%next%ip,iaary(il)%p%next%xb,iaary(il)%p%next%yb,iaary(il)%p%next%s
    !!deb print '(1x,"IB prev", 1024(e16.8,1x))', ibary(il)%p%prev%x, ibary(il)%p%prev%y, ibary(il)%p%prev%s, mkn(1,2,il)
    !!deb print '(1x,"IB     ", 1024(e16.8,1x))', ibary(il)%p%x, ibary(il)%p%y, ibary(il)%p%s, mkn(1,1,il)
    !!deb print '(1x,"IB next", 1024(e16.8,1x))', ibary(il)%p%next%x, ibary(il)%p%next%y, ibary(il)%p%next%s, mkn(1,nring,il)
    !!deb print '(1x,1024(i3,1x))', ibary(il)%p%prev%cor(1,:)
    !!deb print '(1x,1024(i3,1x))', ibary(il)%p%cor(1,:)
    !!deb print '(1x,1024(i3,1x))', ibary(il)%p%next%cor(1,:)
!deb    print *, 'Done for ', il, ' list '
  enddo
!
  return
end subroutine getAugVarslist
!
!
!========================================================================
! find the point of intersection of cell and grid line
!========================================================================
subroutine intercept(xt, yt, nx, ny, s, sz, slow, supp, xl, xu, mk, il, isel, info)
  implicit none
! we have parameter between [slow,supp] and the (xt,yt) is used to find
! intercept: 
! isel = 0, use xt(s) to find parameter slow <= s <=supp and yt(s)
! isel = 1, use yt(s) to find parameter slow <= s <=supp and xt(s)
! here we assume (x(slow),y(slow)) and (x(supp),y(supp)) are on the 
! two sides of (xt(s),yt(s))
  double precision :: xt, yt, nx, ny, xl, xu, s, sz, slow, supp, mk(7,nring,nent-2)
  integer :: il, isel, info, is, js, kid
!
  double precision :: sl, su, sc, sd, xc, yc, eps, vc, tl, tu, tc
  double precision :: pl, pu, pc, tp1, tp2, tmp, ss0, ss1, dth
  integer :: i, imx, ic
!
  imx = 60
  eps = 1.d-8*hg
  ic = isel
  kid = mod(info+nring-2,nring);
  dth = two*cpi/dble(nring);
  ss0 = dth*dble(kid);
!
  select case (isel)
  case (0)
    vc = xt; ! target (middle point)
  case (1)
    vc = yt;
  end select
  pl = xl; pu = xu
!
  sl = slow; su = supp; sd = abs(supp - slow); tl = pl - vc; tu = pu - vc
  if (sl > su) then
    print *, 'wrong calling options stop here'
    stop
  endif
  if (tl*tu > 0) then
    print *, 'wrong searching data, check!'
    print '(1x,1024(e12.5,1x))', slow, s, supp, xl, xu, tl, tu
    stop
  endif
  i = 0
  do while  ( i <= imx .and. sd > eps)
    info = 0
    sc = half*(sl+su); pc = paraval(sc, nring, mk,il,ic,0,info);
    if (info .eq. -1) then
      print '(1x,"check here now ", i5, 1024(e16.8,1x))', info, pc, sd
      print '(1x, 1024(e14.7,1x))', xt, pl, pu, xl, xu, slow, supp, sl, sc, su, tl, tu, tc
      print '(1x, 1024(i5,1x))', i, il, isel, info
    endif
    tc = pc - vc
    if (tc*tl > zero .and. tc*tu < zero) then
      sl = sc; tl = tc; sd = abs(su - sl)
    elseif (tc*tu > zero .and. tc*tl <zero) then
      su = sc; tu = tc; sd = abs(su - sl)
    else
      if (abs(tc) < 1d-15) then
        exit
      endif
      print *, 'cases wrong!'
      print '(1x,1024(e16.8,1x))', xt, pl, pc, pu, xl, xu, slow, supp, sl, slow, su, supp
      print '(1x,1024(e16.8,1x))', dble(i), tl, tc, tu, sl, sc, su, sd
      stop
    endif
    i = i + 1
  enddo
  !!prt print *, 'loop= ', i, 'sd = ', sd
  sz = sc 
!123456789212345678931234567894123456789512345678961234567897123456789
  s  = (sc-slow)/(supp-slow)*dth+ss0
!123456789212345678931234567894123456789512345678961234567897123456789
  xc = paraval(sc, nring, mk, il, ic, 0,info)
  tc = paraval(sc, nring, mk, il, 1-ic, 0,info)
  tp1= paraval(sc, nring, mk, il, 0, 1, info)
  tp2= paraval(sc, nring, mk, il, 1, 1, info)
  tmp= sqrt(tp1*tp1+tp2*tp2)
  nx =-tp2/tmp
  ny = tp1/tmp

  if (ic.eq.0) then
    yt=tc
  else
    xt=tc
  endif
  !!if (abs(xc-tc)>eps) then
  !!  info =-1
  !!  print *, 'This is not the point we are looking for!'
  !!  stop
  !!endif
  if (info .eq. -1) then
    print *, 'check here now part 2'
    print '(1x, 1024(e16.8,1x))', xt, pl, pu, slow, supp, sl, su, tl, tu, tc
    print '(1x, 1024(i5,1x))', il, isel, info
  endif
!
  return
end subroutine intercept
!
!
subroutine cc_errmsg(info, iblst, cp_ss)
  implicit none
  integer :: info
  type(ibpt), dimension(:) :: iblst(nent-2)
  type(cc_bp), pointer :: curr, cb_bc
  type(cc_augvar), pointer :: cp_lp, cp_ss
!
  integer :: i, il, imx
  double precision :: tmp
!
  if (info < 0) then
    imx = 20
    do il =1, nib
      curr => iblst(il)%p
      do i = 1, nring
        print '(1x,1024(e16.8,1x))', curr%x, curr%y
        curr => curr%prev
      enddo
    enddo
    print *, 'finish IB'

    i = 1
    do while ( associated(cp_ss) .and. i .le. imx) 
      print '(1x,1024(e16.8,1x))', cp_ss%xb, cp_ss%yb
      i = i + 1
      cp_ss => cp_ss%prev
    enddo
    call abort
    
  endif

!
  return
end subroutine cc_errmsg
!
!========================================================================
! save stencil at grid crossing
!========================================================================
subroutine savegridcross(idf,xt,yt, nx, ny, s, sl, dxs, ik, jk, ip, jp, cp_ss, info)
  implicit none
  ! save grid crossing and the stencil involved on both sides (+,-)
  ! of (xb,yb). on input, we also know unit outward normal (nx,ny) and
  ! the index of the two grid point next to (xb,yb), we will use + & -
  ! side to label the two points and store their grid index in data
  ! node
  ! on output 
  ! on input: (xt,yt) grid crossing location
  ! on input: (nx,ny) unit normal at (xt,yt) (point from + to - side)
  ! on input: s is the arclength parameter of (xt,yt)
  ! on input: sl is the material parameter 
  ! on input: (ik,jk) is one of the two grid next to (xt,yt)
  !
  integer, intent(in) :: idf(-1:nxb+1,-1:nyb+1)
  double precision :: xt, yt, nx, ny, s, sl, slow, supp, dxs
  integer :: ik, jk, ip, jp, info
  type(cc_augvar), pointer :: cp_lp, cp_ss, cp_tt
!
  double precision :: x0, y0, x1, y1, xc,yc, delm, delp, st_m(4), st_p(4)
  double precision :: p0(2), p1(2), pt(2)
  integer :: i,j, il(2), ir(2), it, jt, is(2), kst(4,2)
!
  it = idf(ik,jk); 
  if (it .eq. -1) then ! il on "-" side; ir on "+" side
    il(1) =ik;    il(2) = jk
    ir(1) =ik+ip; ir(2) = jk+jp
  else if (it.eq. 1) then
    ir(1) =ik;    ir(2) = jk
    il(1) =ik+ip; il(2) = jk+jp
  else
    print *, 'there should be no such option!'
    print '(1x,1024(i3,1x))', ik, jk, ip, jp, it, il, ir
    stop
  endif
  it = idf(il(1),il(2))
  jt = idf(ir(1),ir(2))
  if (it*jt .eq. 1) then
    print *, 'Two ends should be on different sides of boundary! check it'
    print '(1x,1024(i3,1x))', il, ir, it, jt, idf(il(1),il(2)), idf(ir(1),ir(2)), ik, jk, ip, jp
  x0 = xmin+(dble(il(1))-xshift)*hg; ! idf(x0,y0) =-1
  y0 = ymin+(dble(il(2))-yshift)*hg
  x1 = xmin+(dble(ir(1))-xshift)*hg; ! idf(x1,y1) =1
  y1 = ymin+(dble(ir(2))-yshift)*hg
    print '(1x,1024(e13.6,1x))', xt, yt, x0, y0, x1, y1, nx, ny, s
    info = -3
    goto  987
!!    stop
  endif

  x0 = xmin+(dble(il(1))-xshift)*hg; ! idf(x0,y0) =-1
  y0 = ymin+(dble(il(2))-yshift)*hg
  x1 = xmin+(dble(ir(1))-xshift)*hg; ! idf(x1,y1) =1
  y1 = ymin+(dble(ir(2))-yshift)*hg
  if ((y1-yt)*(yt-y0)+(x1-xt)*(xt-x0) <0) then
    print '(1x, 1024(e16.8,1x))', y1, yt, y0, x1, xt, x0, s
    print '(1x, 1024(i5,1x))', ik, jk, ip, jp, ir(:)
    print *, 'wrong grid crossing point!'
    stop
  endif
  !
  delm = abs((xt - x0) + (yt-y0)) ! distance from "-" neighbor to crossing
  delp = abs((x1 - xt) + (y1-yt)) ! distance from "+" neighbor to crossing
  delm = delm/hg;           delp = delp/hg
  !
  cp_ss%im = il; cp_ss%ip = ir
  cp_ss%s  = s;  cp_ss%dxs= dxs; cp_ss%sl = sl
  cp_ss%xb = xt; cp_ss%yb = yt
!!seemsNoUse  cp_ss%cm = 1.; cp_ss%cp = 1.
  cp_ss%nx = nx; cp_ss%ny = ny
  cp_ss%delm = delm; cp_ss%delp = delp;
!
  is = cp_ss%ip - cp_ss%im ! is> 0 means (xt,yt) is on right/above im
  j = 3; kst = 0
  do i = 1, 3
    kst(i,:) = cp_ss%im(:) - is(:)*(i-1) ! kst(1,:) = il
    if (idf(kst(i,1),kst(i,2)) .ne. -1) then
      j = i-1
      exit
    endif
  enddo
  cp_ss%iml= j ! # of grids can be used for extrapolation in irregular grid 
               ! and for directional derivative discretization along grid line
  call stencil_ali_grd(is,j, delm, st_m,info) ! for extrapolation along grid line
  kst(4,1) = j
  cp_ss%ilocm = kst; cp_ss%st_m = st_m;
  !!cp_ss%st_m = st_m ! note st_m only store extrapolation coefficient
!
  info = 0
  !!old version call setgradient(il,ir,xt,yt,nx,ny,delm,delp,ost,pst,-1,j, info)
  call setgradient(idf,cp_ss,-1, info)
  !
  if (info .eq. -1) then
    print *, 'calling directional derivative wrong on - side'
  endif
!
  is = cp_ss%im - cp_ss%ip ! not st_p only store extrapolation coefficient
  j = 3; kst = 0
  do i = 1, 3
    kst(i,:) = cp_ss%ip - is*(i-1)
    if (idf(kst(i,1),kst(i,2)) .ne. 1) then
      j = i-1
      exit
    endif
  enddo
  cp_ss%ipl= j
  call stencil_ali_grd(is,j, delp, st_p,info)
  kst(4,1) = j
  cp_ss%ilocp = kst; cp_ss%st_p = st_p
!
  info = 0
  !!call setgradient(il,ir,xt,yt,nx,ny,delm,delp,ost,pst,1,j, info)
  call setgradient(idf,cp_ss,1, info)
  if (info .eq. -1) then
    print *, 'calling directional derivative wrong on + side'
    print '(1x,1024(i3,1x))', idf(il(1),il(2)), idf(ir(1),ir(2)),  cp_ss%ipl, info
  endif
!
!=======================================================================
! check the discretization of BC at grid crossing are done in 
! savegradient. on each side, we will have two directions: grid aligning
! %c_ali_grid; and off grid %c_off_grid. the actual grid points involved
! should be in %i_ali_grid, and %i_off_grid
!!    double precision :: cm, cp, xb, yb, delm, delp, s, nx, ny
!!    double precision :: st_m(4), st_p(4), locm(6,2), locp(6,2)
!!    integer :: ilocm(4,2),ilocp(4,2)
!!    integer :: im(2), ip(2), iml, ipl!index of - & + side of (xb,yb)
!!    integer :: i_ali_grd_m(4,2), i_off_grd_m(4,2)
!!    integer :: i_ali_grd_p(4,2), i_off_grd_p(4,2)
!!    double precision :: c_ali_grd_m(4), c_off_grd_m(4)
!!    double precision :: c_ali_grd_p(4), c_off_grd_p(4)
!=======================================================================
987 continue
!
  return
end subroutine savegridcross
!
!
subroutine setgradient(idf,cp_ss,isel, info)
  implicit none
! on input,  isel = 1 need to find grids on (ip,jp) side,
!            isel =-1 need to find grids on (im,jm) side
! so that all grid points needed to discretize normal at (xp,yp)
! on output i_ali_grd_?(4,1) give index of grids used along grid, and 
!           i_off_grd_?(4,2) give index of grids used on off grid direction
! the coefficients are stored in c_ali_grd_?(4) and c_ali_grd_?(4)
! isel determ _? to be _m (isel = -1) or _p (isel = 1)
! note cp_ss%ilocm & cp_ss%ilocp should have same content as in 
! cp_ss%i_ali_grd_[m|p], and we may reduce the redundancy in future,
! by removing iloc[m|p], for debugging, p_ali_grd_[m|p] & p_off_grd_[m|p]
! are introduced to check points involved in discretizing BCs at (xp,yp)^[+|-]
!
  integer :: idf(-1:nx+1,-1:ny+1)
  type(cc_augvar), pointer :: cp_lp, cp_ss
  integer :: im(2),ip(2),info, isel, cpl, opl, gpl, jcase, ost(4,2), iloc(4,2)
  double precision :: xp,yp, delm, delp, nx, ny, mv(2), pst(3,2), qst(3,2), igrd
!
  double precision :: tmp, mx, my, delta, c_ali_grd(4), c_off_grd(4)
  integer :: ilocm(4,2), ilocp(4,2)!, i_ali_grd(4,2), i_off_grd(4,2)
  !!double precision :: p_ali_grd(4,2), p_off_grd(4,2)
  integer :: i, j, i0(2), i1(2), i2(2), j2(2), is(2)
  double precision :: smt(2,2), ab(2), ar(2), ugrd(2), nugrd, uoff(2), nuoff, determ
  type(cc_bp), pointer :: curr
  double precision :: x0(2), x1(2), x2(2), xgx(2), rh, sten(4)
!=======================================================================
! read in data
  xp = cp_ss%xb; yp = cp_ss%yb; nx = cp_ss%nx; ny = cp_ss%ny
  im = cp_ss%im; ip = cp_ss%ip; delm = cp_ss%delm; delp = cp_ss%delp
  if (isel .eq. -1) then
    cpl= cp_ss%iml; delta = delm; iloc=cp_ss%ilocm; sten = cp_ss%st_m
  else
    cpl= cp_ss%ipl; delta = delp; iloc=cp_ss%ilocp; sten = cp_ss%st_p
  endif
  do i = 1,cpl
    if (chker1(iloc(i,1),iloc(i,2)) .ne. isel) then
    !!if (idf(iloc(i,1),iloc(i,2)) .ne. isel) then
      print *, 'Wrong case here, need to check'
! debug      
      print '(1x,1024(i3,1x))', isel
      print '(1x,1024(i3,1x))', (iloc(j,:),j=1,3), (idf(iloc(j,1),iloc(j,2)),j=1,3), iloc(4,1)
      stop
    endif
  enddo
  xgx(1) = xp; xgx(2) = yp; rh = 1./hg
!=======================================================================
!
  if (isel .eq. 1) then ! need to find grids on (ip,jp) side
    i0 = im; i1= ip; mx = -nx; my = -ny ! get the normal pointing into + side
  else if (isel .eq. -1) then
    i0 = ip; i1= im; mx =  nx; my =  ny! get the normal pointing into - side
  else
    print *, 'there should be no such option instead of just 1 or -1'
    print *, 'isel = ', isel
    stop
  endif
  !!mv(1) =-mx; mv(2) =-my ! the unit normal is flipped here 
  mv(1) = nx; mv(2) = ny ! the unit normal are the same
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! Finding all the grid points that will be used for off grid direction
!=======================================================================
  info = -1; jcase = 0; ost = 0
  is(1) = int(dsign(1.d0,dble(mx))); is(2) = int(dsign(1.d0,dble(my)))
  i2 = i0+is ! first point
  if (idf(i2(1),i2(2)) .eq. isel) then
  !!if (chker0(i2(1),i2(2)) .eq. isel) then
    ost(1,:) = i2 ! first point
    info = 1
  else 
    print '(1x,"boundary is not resolved properly, only 1 off grid point!", 2(i5,1x))', i0
    ! this case happens if 3 grid cross have 1 common edge. the grid cross 
    ! that is in middle will have to use 
    !!print '(1x,1024(i3,1x))', i0, i1, is, i2, (ost(j,:),j=1,2),&
    !! (idf(ost(j,1),ost(j,2)),j=1,2),idf(i2(1),i2(2)), isel, cpl, info
    !!print '(1x,1024(i3,1x))', i0, i1, is, i2, idf(i0(1),i0(2)), idf(i1(1),i1(2)), idf(i2(1),i2(2)), isel, cpl, info
    !!deb do i = 1, 2
    !!deb   x0(i) =  xmin+h*(dble(i0(i))-xshift)
    !!deb   x1(i) =  xmin+h*(dble(i1(i))-xshift)
    !!deb   x2(i) =  xmin+h*(dble(i2(i))-xshift)
    !!deb enddo
    !!deb print '(1x,1024(e16.8,1x))', x0, xp,yp, x1, x2, mx, my, nx, ny, delm, delp, h
!
    i2 = i1+is
    if (idf(i2(1),i2(2)) .eq. isel) then
      ost(1,:) = i2
      info = 1; jcase = 1 ! special case
      i2 = i1 + 2*is
      if (idf(i2(1),i2(2)) .eq. isel) then
        ost(2,:) = i2
        info = 2
      else
        print *, 'double precisionly can not get high order'
      endif

    else
      print *, 'out of idea, has to stop', isel
      print *, is, i1, i2
      stop
    endif
  endif
  i2 = i0+2*is;
  if (idf(i2(1),i2(2)) .eq. isel) then
    ost(2,:) = i2
    info = 2
    j2 = i2+i0-i1
    if (idf(j2(1),j2(2)) .eq. isel) then
      ost(3,:) = j2
      info = 3
    else
      info = 1
      !!prt print *, 'can only has 1st order approximation'
      !!prt print '(1x,"1st order A", 1024(i3,1x))', i0, i1, i2, is, idf(i2,j2), isel
      !!prt print *, nx, ny, xp, yp
    endif
  else 
    info = 1
!info    print '(1x,"1st order extrapolation ", 10(e14.6,1x))', cp_ss%xb, cp_ss%yb
!    print *, 'can only has 1st order approximation'
!    print '(1x,1024(i3,1x))', i0, i1, is, i2, j2, idf(i2,j2), isel
  endif
  ! set order of interpolation used for directional derivative in opl
  !!opl = 1
  opl = info
  ost(4,1)=opl; ost(4,2)=info! this give # of "info" grid points for off grid pts
  !!ost(4,1)=cpl; ost(4,2)=1! use only 1 grid point
!
  pst = 0.
  !!do i = 1, info ! get locations of all points found, some of them may not be used
  do i = 1, opl ! get locations of all points found, some of them may not be used
    pst(i,1) = xmin+hg*(dble(ost(i,1))-xshift)
    pst(i,2) = ymin+hg*(dble(ost(i,2))-yshift) 
  enddo
  do i = 1, cpl ! get locations of all points found, some of them may not be used
    qst(i,1) = xmin+hg*(dble(iloc(i,1))-xshift)
    qst(i,2) = ymin+hg*(dble(iloc(i,2))-yshift) 
  enddo
!=======================================================================
! need to compute stencil along/off grid line for directional derivative
!
  uoff = 0.; ugrd = 0.
  !uoff(1) = pst(1,1) - xp; uoff(2) = pst(1,2) - yp
  uoff = pst(1,:) - xgx ! first point in off grid stencil
  if (cpl .eq. 3 .or. cpl .eq. 2) then ! at least 1 grid pt along grid line
    !ugrd(1) = qst(2,1) - xp; ugrd(2) = qst(2,2) - yp
    ugrd = qst(2,:) - xgx
  else if(cpl .eq. 1) then ! only irregular grid pt along grid line available
    !ugrd(1) = qst(1,1) - xp; ugrd(2) = qst(1,2) - yp
    ugrd = qst(1,:) - xgx
  else
    print *, 'should not take this value!', cpl
    stop
  endif
  nugrd = sqrt(dot_product(ugrd,ugrd));
  nuoff = sqrt(dot_product(uoff,uoff)); 

  smt(1,1) = dot_product(ugrd, ugrd); smt(1,2) = dot_product(ugrd, uoff);
  smt(2,1) = dot_product(ugrd, uoff); smt(2,2) = dot_product(uoff, uoff); 
  ar(1) = dot_product(mv, ugrd)
  ar(2) = dot_product(mv, uoff)
  determ = smt(1,1)*smt(2,2) - smt(1,2)*smt(2,1)
  !!smt = smt/determ
  !!ar  = ar /determ
  !!determ = 1.
  ab(1) = ( ar(1)*smt(2,2)-ar(2)*smt(1,2))/determ
  ab(2) = (-ar(1)*smt(2,1)+ar(2)*smt(1,1))/determ
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! along grid line direction: note diffusion coeff is not included`
  tmp = ab(1)*nugrd
  if (cpl > 1) then
    gpl = 2; ! gpl=1 1st order everywhere, gpl=cpl whatever available
  else
    gpl = 1
  endif
  gpl = cpl

  iloc(4,1) = gpl !store order choice in grid array
  iloc(4,2) = cpl !store # of grid pts available
  !!if (cpl .eq. 3) then! use all grid points found
  if (gpl .eq. 3) then ! gpl control order ! should use select case
    c_ali_grd(1) =-(3.+2.*delta)/(1.+delta)/(2.+delta)
    c_ali_grd(2) = 0.
    c_ali_grd(3) = (2.+delta)/(1.+delta)
    c_ali_grd(4) =-(1.+delta)/(2.+delta)
    c_ali_grd    = tmp*rh*c_ali_grd
    !old c_ali_grd(1) = sten(1)
    !old c_ali_grd(2) = -1.
    !old c_ali_grd(3) = sten(3)
    !old c_ali_grd(4) = sten(4)
    !old c_ali_grd    =-ab(1)*c_ali_grd*nugrd*rh
  !else if (cpl .eq. 1) then! use all grd pts available
  else if (gpl .eq. 2) then ! gpl control order
    c_ali_grd(1) =-1. ! coefficient at grid crossing
    c_ali_grd(2) = 0.    ! this is the special case
    c_ali_grd(3) = 1.
    c_ali_grd(4) = 0.
    c_ali_grd    = c_ali_grd*ab(1)

    !!c_ali_grd(1) = 2./(1.+delta) ! coefficient at grid crossing
    !!c_ali_grd(2) =-1.    ! this is the special case
    !!c_ali_grd(3) =-(1.-delta)/(1.+delta)
    !!c_ali_grd(4) = 0.
    !!c_ali_grd    =-c_ali_grd*ab(1)*nugrd*rh
    !!ck c_ali_grd(1) =-1. ! coefficient at grid crossing
    !!ck c_ali_grd(2) = 1.-delta    ! this is the special case
    !!ck c_ali_grd(3) = delta
    !!ck c_ali_grd(4) = 0.
    !!c_ali_grd    = c_ali_grd*nugrd*rh*ab(1)
  else if (gpl .eq. 1) then ! gpl control order
    c_ali_grd(1) =-ab(1) ! coefficient at grid crossing
    c_ali_grd(2) = ab(1) ! this is the special case
    c_ali_grd(3) = 0.
    c_ali_grd(4) = 0.
    !!c_ali_grd    = c_ali_grd
  else
    print *, 'we did not prepare for this case', cpl
    stop
  endif
  !!c_ali_grd = c_ali_grd*dif(1)
  c_ali_grd = c_ali_grd
! off grid line direction: note diffusion coeff is not included`
  !!if (info .eq. 3) then ! use all available grid pts found
  tmp = ab(2)*nuoff
  if (opl .eq. 3) then ! opl control order of interpolation
    c_off_grd(1) =-1.5
    c_off_grd(2) = 2.
    c_off_grd(3) =-0.5*delta
    c_off_grd(4) =-0.5*(1.-delta)
    c_off_grd    = c_off_grd*ab(2)
  !!else if (info .eq. 1) then! use all available grid points found
  else if (opl .eq. 1 .or. opl .eq. 2) then! opl control order of interpolation
    c_off_grd(1) =-1. 
    c_off_grd(2) = 1.
    c_off_grd(3) = 0.
    c_off_grd(4) = 0.
    c_off_grd    = c_off_grd*ab(2)
  else
    print *, 'off grid line, we did not prepare for this', opl
    stop
  endif
  !!c_off_grd = c_off_grd*dif(1)
  c_off_grd = c_off_grd
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  !!prt print '(1x,"directions! ", 1024(e16.8,1x))', ugrd, uoff
!=======================================================================
!=======================================================================
! save grid indexing info to node
  if (isel .eq. -1) then
    cp_ss%i_off_grd_m = ost; cp_ss%i_ali_grd_m = iloc;
    cp_ss%p_off_grd_m = pst; cp_ss%p_ali_grd_m = qst
    cp_ss%c_off_grd_m = c_off_grd;
    cp_ss%c_ali_grd_m = c_ali_grd
  elseif (isel .eq. 1) then
    cp_ss%i_off_grd_p = ost; cp_ss%i_ali_grd_p = iloc;
    cp_ss%p_off_grd_p = pst; cp_ss%p_ali_grd_p = qst
    cp_ss%c_off_grd_p = c_off_grd;
    cp_ss%c_ali_grd_p = c_ali_grd
  else
    print *, 'no such option'
    stop
  endif
!=======================================================================
!
  return
end subroutine setgradient
!
subroutine stencil_ali_grd(is,j,delp,st_p,info)
  double precision :: st_p(4), delp
  integer :: is(2), j, info
! THE TWO HAVE NO DIFFERENCE AND SHOULD BE COMBINED 
!!  if (sum(is) > 0) then
    select case (j) 
    case (3)
      st_p(1) = 6./(1.+delp)/(2.+delp)
      st_p(2) = 0.
      st_p(3) =-3.*(1.-delp)/(1.+delp)
      st_p(4) = 2.*(1.-delp)/(2.+delp)
    case (2)
      st_p(1) = 2./(1.+delp)
      st_p(2) = 0.
      st_p(3) =-(1.-delp)/(1.+delp)
      st_p(4) = 0.
    case (1)
      st_p(1) = 1.
      st_p(2) = 0.
      st_p(3) = 0.
      st_p(4) = 0.
    case default
      print *, 'boundary is not resolved B!'
      stop
    end select 
!
  return
end subroutine stencil_ali_grd

!
!123456789212345678931234567894123456789512345678961234567897123456789

!subroutine smoothcrossing(m,ni,il,mkc,iaary,llen,md,pd,info)
subroutine smoothcrossing(m,mkc,iaary,llen,md,pd,info)
! on input mkc(1,1:nring, :) store all arclength parameters @ IBs
  integer :: m(nent-2), llen(nent-2), info
  double precision :: mkc(7,nring,nent-2), mkx(7,nring,nent-2)
  double precision :: md(nring,nent-2), pd(nring,nent-2) ! only md(1:llen,nid) & pd(1:llen,nid) store data
  type(iapt), dimension(:) :: iaary(nent-2)
  type(ibpt), dimension(:) :: iblst(nent-2)
!
  integer :: n, i, j, k, il, ipiv(3), iwork(6), mct
  !double precision :: tmp, mI(3,3), mId(3,3), s(nring+1), mA(3,3), xx(3), pc(3), &
  double precision :: tmp, s(nring+1), mA(3,3), xx(3), pc(3), &
    rcond, mC(3), mR(3), ferr(1), berr(1), AF(3,3), work(24), rx(3), pA(3,3)
  double precision :: s0, s1, s2, tupi, x0, x1, x2, y0, y1, y2
  double precision :: ffi(4*nring+15)!, mdi(n), pdi(n), smd(n), spd(n), sw(n)
  type(cc_augvar), pointer :: cp_lp, cp_ss
  character :: equed
  double precision :: cc_tp(nring)
  integer :: lensav, ier, lenwrk, inc
  double precision :: cc_cm(nring,nent-2), cc_cp(nring,nent-2)
  integer*8 :: plan
  complex*16, dimension(nring/2+1) :: dout
!
  !!if (n .ne. llen(il)) then
  !!  print *, 'inconsistent input!'
  !!  stop
  !!endif
  n = nring
!
  tupi = two*cpi
!  mI = zero; 
!  do i = 1, 3
!    mI(i,i) = one
!  enddo

  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      mkx(1,i,il) = cp_lp%s; 
!!    md(i,il) = cp_lp%cm; 
!!    pd(i,il) = cp_lp%cp;
      cp_lp => cp_lp%next
    enddo
  enddo

  do il = 1, nib
    call Bubble_Sort(llen(il),nib,il,mkx,md(:,il),pd(:,il),info)
    call removedup(llen(il), m(il), il, mkx, md(:,il), pd(:,il), info)
!deb    print *, 'in smoothcrossing', nring, m(il), llen(il)
  enddo

  do il = 1, nib
    do i = 1, nring
      cc_cm(i,il) = linint(mkc(1,i,il), m(il), nring, mkx, il, md(:,il), 0, info)
      cc_cp(i,il) = linint(mkc(1,i,il), m(il), nring, mkx, il, pd(:,il), 0, info)
    enddo
    !!prt print *, ' record '
    !!prt do i = 1, nring
    !!prt   print '(1x,i3,1x,1024(e16.8,1x))', i, mkc(1,i,il), cc_cm(i,il), cc_cp(i,il)
    !!prt enddo
    !!prt print *, ' record done'
    mct = (nring)/16
    cc_tp = cc_cm(:,il)
    call dfftw_plan_dft_r2c_1d(plan,nring,cc_tp,dout,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, cc_tp, dout)
    call dfftw_destroy_plan(plan)
    dout(nring/2+2-mct:nring/2+1) = cmplx(zero,zero)
    call dfftw_plan_dft_c2r_1d(plan,nring,dout,cc_tp,FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan, dout, cc_tp)
    cc_cm(:,il) = cc_tp/dble(nring)
    call dfftw_destroy_plan(plan)

    cc_tp = cc_cp(:,il)
    call dfftw_plan_dft_r2c_1d(plan,nring,cc_tp,dout,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, cc_tp, dout)
    call dfftw_destroy_plan(plan)
    dout(nring/2+2-mct:nring/2+1) = cmplx(zero,zero)
    call dfftw_plan_dft_c2r_1d(plan,nring,dout,cc_tp,FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan, dout, cc_tp)
    cc_cp(:,il) = cc_tp/dble(nring)
    call dfftw_destroy_plan(plan)
  enddo
!=======================================================================
! done smoothing
!=======================================================================
  do il = 1, nib
  mA(:,1) = one
  do k = 1, nring
  !!do k = 1, m(il)
    if (k .eq. 1) then
      s0 = mkc(1,n,il); s1 = mkc(1,k,il)+tupi; s2 = mkc(1,k+1,il)+tupi
      x0 = cc_cm(n,il); x1 = cc_cm(k,il); x2 = cc_cm(k+1,il)
      y0 = cc_cp(n,il); y1 = cc_cp(k,il); y2 = cc_cp(k+1,il)
      !!x0 = md(n); x1 = md(k); x2 = md(k+1)
      !!y0 = pd(n); y1 = pd(k); y2 = pd(k+1)
    elseif (k .eq. nring) then
    !!elseif (k .eq. m(il)) then
      !!s0 = mkc(1,k-1,il); s1 = mkc(1,k,il); s2 = 0.
      s0 = mkc(1,k-1,il); s1 = mkc(1,k,il); s2 = tupi
      x0 = cc_cm(k-1,il); x1 = cc_cm(k,il); x2 = cc_cm(1,il);
      y0 = cc_cp(k-1,il); y1 = cc_cp(k,il); y2 = cc_cp(1,il);
      !!x0 = md(k-1); x1 = md(k); x2 = md(1);
      !!y0 = pd(k-1); y1 = pd(k); y2 = pd(1);
    else
      s0 = mkc(1,k-1,il); s1 = mkc(1,k,il); s2 = mkc(1,k+1,il)
      x0 = cc_cm(k-1,il); x1 = cc_cm(k,il); x2 = cc_cm(k+1,il)
      y0 = cc_cp(k-1,il); y1 = cc_cp(k,il); y2 = cc_cp(k+1,il)
      !!x0 = md(k-1); x1 = md(k); x2 = md(k+1)
      !!y0 = pd(k-1); y1 = pd(k); y2 = pd(k+1)
    endif
    mA(1,2) = s0; mA(1,3) = s0**2.d0;
    mA(2,2) = s1; mA(2,3) = s1**2.d0;
    mA(3,2) = s2; mA(3,3) = s2**2.d0;
    rx(1) = x0; rx(2) = x1; rx(3) = x2
    pc = rx
    pA = mA
    CALL DGESVX('Equilibration','No transpose',3,1,pA,3,AF,3,ipiv,   &
           equed,mR,mC,pc,3,xx,3,RCOND,FERR,BERR,WORK,IWORK,INFO)
    if (info .ne. 0) then
      print *, k, n,m, llen(1), nring
      print '(1x,"case x", 12(e16.8,1x))', rcond, mkc(1,k-1,il),mkc(1,k,il), mkc(1,k+1,il)
      print '(1x,1024(e16.8,1x))', s0, s1, s2
      print '(1x,1024(e16.8,1x))', x0, x1, x2
      print '(1x,1024(e16.8,1x))', rx
      print '(1x,1024(e16.8,1x))', pc
      print *, ''
      do i = 1, 3
        print '(1x,1024(e16.8,1x))', mA(:,i)
      enddo
      print *, 'above mA'
      do i = 1, 3
        print '(1x,1024(e16.8,1x))', pA(:,i)
      enddo
      print *, 'above pA'
      print '(1x,1024(e16.8,1x))', xx
      cp_lp => iaary(1)%p
      do i = 1, llen(1)
        print '(1x,5(i3, 1x),1024(e16.9,1x))', i, cp_lp%im, cp_lp%ip, cp_lp%s, cp_lp%xb, cp_lp%yb,&
         mkc(1,i,1), md(i,il), pd(i,il)
        cp_lp =>cp_lp%next
      enddo
      stop
    endif
    mkc(2:4,k,il) = xx

    rx(1) = y0; rx(2) = y1; rx(3) = y2
    pc = rx
    pA = mA
    CALL DGESVX('Equilibration','No transpose',3,1,pA,3,AF,3,ipiv,   &
           equed,mR,mC,pc,3,xx,3,RCOND,FERR,BERR,WORK,IWORK,INFO)
    if (info .ne. 0) then
      print '(1x,"case y",12(e13.6,1x), "case y")', rcond
      stop
    endif
    mkc(5:7,k,il) = xx
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  enddo
  enddo

!
  return
end subroutine smoothcrossing
!
subroutine removedup(n, m, il, mkA, md, pd, info)
  integer :: n, il, m, info
  double precision, intent(in out) :: mkA(7,nring,nib), md(nring), pd(nring)
!
  double precision :: tmp, tp1, tpv(7), tm, tp, eps
  integer :: i, j, ict
!
  !!eps = 1d-1*2.*cpi/dble(nring)
  eps = 0.1d0*two*cpi/dble(nring)
!  do il = 1, nib
    ict = 0
    do i = 2, n-1! not test the last one, may be wrong
      if (abs(mkA(1,i,il)-mkA(1,i-1,il)) < eps) then
        ict = ict + 1
        tpv = mkA(:,i,il); tm = md(i); tp = pd(i)
        mkA(:,i:n-ict,il) = mkA(:,i+1:n-ict+1,il)
        md(i:n-ict) = md(i+1:n-ict+1)
        pd(i:n-ict) = pd(i+1:n-ict+1)
        mkA(:,n-ict+1,il) = tpv; md(n-ict+1) = tm; pd(n-ict+1) = tp;
      endif
    enddo
    m = n - ict
!  enddo
!
  return
end subroutine removedup
!=======================================================================
subroutine Bubble_Sort(n, ni, il, mkA, md, pd, info)
! bubble sort to order rows in mkA, according to the first column of mkA
  integer :: n, ni, il, info
  double precision, intent(in out) :: mkA(7,nring,ni), md(nring), pd(nring)
  !!double precision, intent(in out), dimension(:) :: a
  double precision :: tpv(7), tm, tp
  integer :: i, j
  logical :: swapped 
!
  !!!do j = size(a)-1, 1, -1
  do j = n-1, 1, -1
    swapped = .false.
    do i = 1, j
      if (mkA(1,i,il) > mkA(1,i+1,il)) then
        !!temp = a(i)
        tpv(:) = mkA(:,i,il)
        tm = md(i); tp = pd(i)
        !!a(i) = a(i+1)
        mkA(:,i,il) = mkA(:,i+1,il)
        md(i) = md(i+1); pd(i) = pd(i+1)
        !!a(i+1) = temp
        mkA(:,i+1,il)=tpv(:)
        md(i+1) = tm; pd(i+1) = tp
        swapped = .TRUE.
      end if
    end do
    if (.not. swapped) exit
  end do
!
  RETURN
end subroutine bubble_sort
!=======================================================================
!
double precision function linint(s, lda, n, mk, il, datv, icase, info)
! on input: n total # of record in mk & datv; only first lda records are used, 
! on output: linint(s) return linear interpolation datv data @ s
! note 1) datv(1:lda) are supposed to be periodic;
!      2) mk(1:m) increases in values
!      3) mk(1:m) are bounded in [0, 2pi]
!      4) assume  0<= s<=2*pi; and mk(1,1,il) > 0
  integer :: lda, n
  double precision :: s, mk(7,n, nent-2), datv(n)
  integer :: il, icase, info
!
  double precision :: ck, tupi, eps, xl, xu, sl, su, tmp
  integer :: i, k, kl, ku, inc, ijump
! 
  tupi = 2.*cpi
!
  eps = 1d-10
  sl = zero; su = zero; xl = zero; xu = zero;
  info = -1
  ! locate k
  ck = mod(mod(s, tupi)+tupi, tupi) ! archlength we will used
!  ck = s0
!PRT  print *, 's0', s0, 'ck', ck
  k = 1
  ijump = 0
  if ( ck > eps ) then
    if (ck < mk(1,k,il)) then
      ck = ck + tupi ! 
      ijump = 1 ! here k = 1
      sl = mk(1,lda,il); su = tupi + mk(1,1,il)
      xl = datv(lda);    xu = datv(1)
      goto 754
    endif
    do while (k < lda .and. ijump .eq. 0)
      if ( ck .ge. mk(1,k,il) .and. ck < mk(1,k+1,il)) then
        ! this does not count case: 2*pi> ck > mk(1, lda, il) 
        ijump = 1 ! 
        sl = mk(1,k,il); su = mk(1,k+1,il)
        xl = datv(k);    xu = datv(k+1)
      else
        k = k + 1 ! ckpt 1
      endif
    enddo
  else ! ckpt 0
    !!ck = ck + tupi
    k = 1
    ijump = 1
  endif
  if (k .eq. lda) then 
    if (ijump .eq. 0) then ! come from ckpt 1
      if ( ck .ge. mk(1,k,il) .and. ck < tupi) then! check if 2*pi> ck > mk(1, lda, il) 
        sl = mk(1,k,il); su = tupi + mk(1,1,il)
        xl = datv(k);    xu = datv(1)
        ijump = 1
      else
        print '(1x,"double precisionly?", 10(e16.8,1x))', ck, s, mk(1,1,il), mk(1,lda,il)
        stop
      endif
    else
      print *, 'this could not happen! check carefully'
      stop
    endif
  endif
  if (k .eq. 1 .and. ijump .eq. 1) then ! this means ck = 0. & from ckpt 0
    !!print *, 'small ck', ck
    ck = tupi + ck ! this is here only if mk(1,1,il)=0
    !!ck = mod(ck, tupi)+tupi ! archlength we will used
    sl = mk(1,lda,il); su = tupi + mk(1,1,il)
    xl = datv(lda);    xu = datv(1)
  endif
  754 continue
!
  tmp = su - sl
  if (abs(tmp) > eps) then
    linint = (xu-xl)/tmp*(ck - sl) + xl
  else
    print *, 'two grid crossing are having same parameter?'
    linint = -1.
    stop
  endif
!
!
  return
end function linint
!
!========================================================================
subroutine releasemem(ibary,iaary,llen)
  type(ibpt), dimension(:) :: ibary(nent-2)
  type(iapt), dimension(:) :: iaary(nent-2)
  integer :: llen(nent-2), isel
!
  integer :: il, i, ierr
  type(cc_bp), pointer :: cb_lp, cb_pp
  type(cc_augvar), pointer :: cp_lp, cp_lq
!
    do il = 1, nib
      cb_lp => ibary(il)%p
      nullify(ibary(il)%p)
      do i = 1, nring
        cb_pp => cb_lp; nullify(cb_lp%next)
        cb_lp => cb_lp%prev
        deallocate(cb_pp)
      enddo
    enddo
    do il = 1, nib
      cp_lp => iaary(il)%p
      nullify(iaary(il)%p)
      i = 0
      do i = 1, llen(il)
        cp_lq => cp_lp;
        if(associated(cp_lp%next)) nullify(cp_lp%next)
        cp_lp => cp_lp%prev
        deallocate(cp_lq, stat=ierr)
!P        print *, ' Status ', i, il, ilp(il), ierr
      enddo
    enddo
!
  return
end subroutine releasemem
!========================================================================
subroutine getGeometry(ibary,iaary, llen,isel, uin, dt)
  type(ibpt), dimension(:), intent(out) :: ibary(nent-2)
  type(iapt), dimension(:), intent(out) :: iaary(nent-2)
  double precision :: dt
  double precision, dimension(:,:,:) :: mk(7,nring,nent-2), nb(nring,2, nent-2)
  double precision, dimension(:,:,:) :: uin(2,nring,nent-2)
  integer :: llen(nent-2), isel
  type(cc_bp), pointer :: curr

  integer :: i,j, im,jm, il
!!
  call getIBlist(ibary,dt)

!========================================================================
! To setup bdry representation by quadratic polynomial, and get unit
! normal at each IB point, mk() & nb() should have same ordering of IB pts
! uin stores dx/dt-u at each IB point
  call BdyQuadParametric(ibary,mk(:,:,1:nent-2),nb, uin)
  ndv(:,:,1:nent-2) = nb(:,:,1:nent-2) ! save to global variable ndv
!!  call BdyCubic(ibary,mk,ccpc,nb,uin,isel)
  if (isel .eq. -1) then ! save to global variable mkn & mko
    mko(:,:,1:nent-2) = mk(:,:,1:nent-2)
    mkn(:,:,1:nent-2) = mk(:,:,1:nent-2)
  else
    mko(:,:,1:nent-2) = mkn(:,:,1:nent-2)
    mkn(:,:,1:nent-2) = mk(:,:,1:nent-2)
  endif
!  call setidx(idx,cip,cbt,ctp,clf,crt,cbl,cbr,ctl,ctr) ! for chemical computing
!========================================================================
! for isel .eq. -1, i.e., the initial step, only chker0 has valid info, so
! idn is set to 0; On exit, chker0=chker1 is the initial tagging & idf=kdf
!========================================================================
! in step 1 after tagging, chker0 has current identifier array, so on exit
! idn has the difference array between first step and initial setup, which
! is to be used for setting RHS in solving chemical.
! while computing idn, chker1 has step 0
! info.  On exit of the first step, chker2 will have the identifier
! array for step 0 and chker1 will have array for step 1 info
!
! for the second step and afterwards, chker 0 have current identifier array,
! so idn will have the difference array for ccpc to use.
! The idk will have the Difference array for ccpo to use. On exit, chker0 will
! store step 2 infor, chker1 will store step 2 info too, since chker0=chker1,
! and chker2 will store step 1 info. In general, after this step, on exit
! of n steps, chker0 will store step n info, chker1 will store n info also,
! and chker2 will store n-1 info. 
!///////////////////////////////////////////////////////////////////////
! for multiple IB objects, need to change Tagging subroutine 
!-----B--1----+----2----+----3----+----4----+----5----+----6----+----7-E

!  call cpu_time(time0)
  call Tagging(ibary,id, idf, oid, kdf, mk, nb, isel)
!  call cpu_time(time1)
!  print *, 'Time used in tagging', time1-time0
!=====================================================

!-----B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
  chker0=idf! set CURRENT checker array to identify pt in/out bdry (1/-1) 
!
  oid = id
  kdf = idf
  ! idn is for checking freshly cleared grid pts
  if (isel .ne. -1)then 
    idn = chker0 - chker1
  else
    idn = 0
  endif
!!  chker2=chker1 ! for level n-2
  chker1=chker0

!!  if (isel .eq. -1) then
!!    chker2=chker0
!!  endif
! set indicator (icor) for ibary
!-----B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
  do il = 1, nib
    curr => ibary(il)%p 
    do i = 1, nring
      do j = 1, 4
        im = curr%cor(j,1); jm = curr%cor(j,2)
        if (id(im,jm) .eq. 1) then
          !id(im,jm) = 1
!          id(im,jm) = il
          curr%icor(j) = 1
        else if (id(im,jm) .eq. -1) then
          curr%icor(j) =-1
        else
          curr%icor(j) = 0
        endif
      enddo
!PPT      print '(1x, 5(i3,1x),2(e13.6,1x))',i,curr%icor,curr%x,curr%y
      curr => curr%prev
    enddo
  enddo
!=======================================================================
  ! The legacy cc_IBbdy path used an ordered list of exact grid/interface
  ! crossings here.  That path is no longer active: the chemical and actin
  ! solvers construct nodewise corrections from side tags and the closest
  ! marker in dualchem/geometry_mod.f90.  Do not build the unused crossing
  ! list, whose one-cell-per-marker-segment assumption rejects otherwise
  ! valid interfaces when adjacent markers map two grid indices apart.
  llen = 0
  do il = 1, nent-2
    nullify(iaary(il)%p)
  enddo
!========================================================================
!
  return
end subroutine getGeometry
!
double precision function mydot(x,y)
  double precision, dimension(-1:nxb+1,-1:nyb+1) :: x, y
!
  integer :: i, j
  double precision :: tsum 
  tsum = 0.d0

  do i = 1, nxb
  do j = 1, nyb
    tsum = tsum + x(i,j)*y(i,j)
  enddo
  enddo
  mydot = tsum
! 
  return
end function mydot
!
!========================================================================
subroutine brs(n,ni,id,mk,cb,s,bp, k,time,isel)
! return the interpolation value of bp(s), with arclength variable provided
! in mk, and function values at corresponding pts in cb, k identifies the seg
! when isel .ne. 0, provide the Dirichlet coeff (u,v)\cdot (nx,ny)
  implicit none
  integer :: n, ni, id, k, isel
  double precision :: bp, s, mk(7,nring,nent-2),cb(nring,nent-2), time
  integer :: i, ijump, info, ipiv(3), iwork(6)
  character :: equed
!
  double precision :: xx(3), mA(3,3), yy(3), mR(3), mC(3), pA(3,3), AF(3,3)
  double precision :: ferr(1), berr(1), work(12)
!!  double precision :: ferr, berr, work(12)
  double precision :: rcond, tupi, s0, s1, s2, tp1, tp2, tmp, x1, y1, vx, vy, x, y
!
  if (isel .ne. 1) then
    tupi = 2.*cpi
    mA(:,1) = 1.
    if (k .eq. 1) then
      s0 = mk(1,n,id) ; 
      s1 = mk(1,1,id) + tupi; 
      s2 = mk(1,2,id) + tupi;
    else if (k .eq. n) then
      s0 = mk(1,k-1,id); s1 = mk(1,k,id); s2 = zero
    else
      s0 = mk(1,k-1,id); s1 = mk(1,k,id); s2 = mk(1,k+1,id)
    endif
    mA(1,2) = s0; mA(1,3) = s0*s0;
    mA(2,2) = s1; mA(2,3) = s1*s1;
    mA(3,2) = s2; mA(3,3) = s2*s2;
!
    xx = zero; yy = zero;
    if ( k .eq. 1) then
      xx(1) = cb(n,id); xx(2) = cb(k,id); xx(3) = cb(k+1,id)
    else if ( k .eq. n) then
      xx(1) = cb(k-1,id); xx(2) = cb(k,id); xx(3) = cb(1,id)
    else
      xx(1) = cb(k-1,id); xx(2) = cb(k,id); xx(3) = cb(k+1,id)
    endif
!
    pA = mA; 
    CALL DGESVX('Equilibration','No transpose',3,1,pA,3,AF,3,ipiv,      &
           equed,mR,mC,xx,3,yy,3,RCOND,FERR,BERR,WORK,IWORK,INFO)
    if (info .ne. 0) then
      print '(1x,"error quit 1", 12(e13.6,1x))', rcond
      stop
    endif 
!
    bp = yy(1) + yy(2)*s + yy(3)*s*s ! time level n-1
  else
    print *, 'not doing this any more'
    stop
  !!  tp1 = 2.d0*mk(4,k,id)*s + mk(3,k,id); 
  !!  tp2 = 2.d0*mk(7,k,id)*s + mk(6,k,id)
  !!  tmp = sqrt(tp1*tp1+tp2*tp2)
  !!  x1  = -tp2/tmp;  y1 =  tp1/tmp; ! unit normal at (xm,ym)
  !!  call crs(nring,nib,id,mk,s,x,y,k) ! find bdry pt for (i,j)
  !!  call cc_vel(x, y, vx, vy, time)
  !!  bp = vx*x1+vy*y1
  endif
!
  return
end subroutine brs
!
! 
subroutine sub2ind(m, i, j, l) 
implicit none 
integer :: m, i, j, l 
 
l=(j-1)*m+i 
 
end subroutine sub2ind
!
subroutine setidx(idx,cip,cbt,ctp,clf,crt,cbl,cbr,ctl,ctr)
! ordering cell centers and marking interior centers, domain edges/corners
  integer, dimension(1:nx*ny) :: idx ! all centers ordered
  integer, dimension(1:(nx-2)*(ny-2)) :: cip ! all interior centers
  integer :: ctl, ctr, cbl, cbr !corners @ cell centers
  integer, dimension(ny-2) :: clf, crt !left, right
  integer, dimension(nx-2) :: cbt, ctp !bottom, top
!
  integer :: i, j, l
!
  l = nx*ny
  do i=1,l
    idx(i)=i
  enddo
!
  do i = 1, ny-2
    cip((i-1)*(nx-2)+1:i*(nx-2))=idx(i*nx+1:(i+1)*nx-1)
  enddo
!!  clf = (/ (i,i= nx+1,l-nx,nx) /)
  do j=1,ny-2
    clf(j)=idx(j*nx+1)
  enddo
!!  crt = (/ (i,i= 2*nx,l-nx,nx) /)
  do j=1,ny-2
    crt(j)=idx((j+1)*nx)
  enddo
!!  cbt = (/ (i, i=2,nx-1) /)
  cbt(1:nx-2)=idx(2:nx-1)
!!  ctp = (/ (i, i=l-nx+1,l-1) /)
  ctp(1:nx-2)=idx(l-nx+1:l-1)
! 4 corners at cell center 
  ctl = l-nx+1
  ctr = l
  cbl = 1
  cbr = nx
!
  return
end subroutine setidx
!
subroutine extend(iaary,llen,mk,udat)
  type(iapt), dimension(:) :: iaary(nent-2)
  !!double precision, dimension(:,:,:) :: uin(2,nring,nent-2)
  double precision, intent (in)  :: mk(7,nring,nent-2)
  integer, intent(in) :: llen(nent-2)
  double precision, dimension(-1:nx+1,-1:ny+1) :: udat
!
  type(cc_augvar), pointer :: cp_lp
  integer :: il, i, j, cpl, info, im(2), ip(2), iloc(4,2)
  double precision :: s0, sbv, tmp
!
  do il = 1, nib
    cp_lp => iaary(il)%p
    do i = 1, llen(il)
      s0 = cp_lp%s; cpl = cp_lp%ipl; iloc = cp_lp%ilocp; im = cp_lp%im
      sbv= paraval(s0, nring, mk, il, 1, 0, info)
!      
      tmp = zero
      do j = 1, cpl
        tmp = tmp + udat(iloc(j,1),iloc(j,2))*cp_lp%st_p(j+1)
      enddo
      tmp = tmp + cp_lp%st_p(1)*sbv
!
      udat(im(1),im(2)) = tmp ! extend network to outside
!
      cp_lp => cp_lp%next
    enddo
  enddo
!
  return
end subroutine extend
!

end module geometry
