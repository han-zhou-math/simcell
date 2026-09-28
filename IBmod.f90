!> module for IB geometry: i) reading/initializing IB pts (InitLocCell);
!>            ii) spread Lagrangian force to Eulerian grid (newSpread);
!> .          iii) interpolate Eulerian variable to IB pts newInterpS
!>
!> NOTE: the implicit dX/dt solve (formerly petmove/linMatMul here, a simpler
!> backward-Euler/single-solve scheme) now lives in fsisolve.f90's AdvanceFSI/
!> zMatMul (midpoint-rule, Picard-iterated, full fluid+elastic coupling).
!>
!>  1) In setting up IB entity (IB object), we assume material coordinates are in [0,2pi]
!>  2) In (xpt,ypt), IB points on single IB object are stored continuously which takes elements, and
!>       are arranged according to the order of IB objects (and fixed after initialization). It is
!>       assumed all IB objects have same # of IB pts
!>  3) The computed force fs(1:mcoor) and scaled force cfs(1:mcoor) store forces @ each IB pts
!>     continuously: (f1_x,f1_y,f2_x,f2_y,.....)
!>  4) All IB objects are assumed closed, i.e., #1 IB pt is neighbor to #2 and #nring IB pts
!> 
!> @author Lingxing Yao
!> @version 0.1
!>  
!> 
!>
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscsnes.h"
module IBmod
  use, intrinsic :: iso_c_binding
    use petsc
    use petscsys
    use petscksp
    use petscsnes
    use petscvec
  use parameters
  use IBforce
!  use geometry, only : paraval
  use myfft

  implicit none
!
private
  double precision, parameter :: h = hg, rh=one/hg
  ! private variables used for updating IB pts (xpt,ypt) in solving dX/dt
  double precision :: vc_nn(mcoor) ! normal direction
  double precision :: vc_sm(mcoor)
  double precision :: vc_cu(mcoor)
  double precision :: vc_nt(mcoor)
  double precision :: vc_us(mcoor)
  double precision :: vc_xo(mcoor)
  double precision :: vc_xs(mcoor)
  double precision :: vc_xt(mcoor)
  ! components stored in above vectors are arranged differently from the force vector fs or cfs

public :: InitLocCell, newInterpS, newSpread, getNormal, curvelen, evalmpoly, &
        & curvePoly, hadamard, normalproject, pusher

contains
!=======================================================================
subroutine curvelen(xb,yb, crvlen)
  double precision, dimension(mpts), intent(in) :: xb,yb
  double precision, intent(out) :: crvlen
!
  integer :: i
  double precision :: tmp
!
  crvlen = zero
  do i = 1, nring-1
    tmp = sqrt((xb(i+1)-xb(i))*(xb(i+1)-xb(i)) + (yb(i+1)-yb(i))*(yb(i+1)-yb(i))) !distance
    crvlen = crvlen + tmp
  enddo
  crvlen = crvlen + sqrt((xb(1)-xb(nring))*(xb(1)-xb(nring))+(yb(1)-yb(nring))*(yb(1)-yb(nring)))
!
  return
end subroutine curvelen
subroutine getNormal(xb,yb, nvx,nvy, mpoly)
!  normal vector (nvx,nvy), 
!  paramatrization mpoly(2:4,i) for x; mpoly(5:7,i) for y at IB index i
! mpoly has parameter in mpoly(1,:) in [0,2pi]
!  crvlen return length of IB polygon (not smoothed)
  double precision, dimension(mpts), intent(in) :: xb,yb
  double precision, dimension(mpts), intent(out) :: nvx,nvy
  double precision, intent(out) :: mpoly(7,nring)
!
  integer :: i, info, ipiv(3), iwork(6)
  double precision :: tmp, spar(mpts), rcond, ferr(1), berr(1), work(24)
  double precision, dimension(3,3) :: mA, AF, pA
  double precision, dimension(3) :: mC, mR, xx, xv, yv, sv, pc
  character :: equed
  !
  spar(1)=zero
  do i = 1, nring-1
    tmp = sqrt((xb(i+1)-xb(i))*(xb(i+1)-xb(i)) + (yb(i+1)-yb(i))*(yb(i+1)-yb(i))) !distance
    spar(i+1) = spar(i) + tmp
  enddo
  tmp = spar(nring) + sqrt((xb(1)-xb(nring))*(xb(1)-xb(nring))+(yb(1)-yb(nring))*(yb(1)-yb(nring)))
  spar = spar*2.d0*cpi/tmp

  mpoly(1,:)=spar(:)
  mA(:,1)=one
  do i = 1, nring
    if (i.eq. 1) then
      sv(1)=spar(nring); sv(2)=two*cpi; sv(3)=spar(2)+two*cpi
      xv(3)=xb(nring);   xv(2)=xb(1);   xv(1)=xb(2);
      yv(3)=yb(nring);   yv(2)=yb(1);   yv(1)=yb(2);
    else if (i.eq. nring) then
      sv(1)=spar(i-1); sv(2)=spar(i); sv(3)=two*cpi
      xv(3)=xb(i-1);   xv(2)=xb(i);   xv(1)=xb(1);
      yv(3)=yb(i-1);   yv(2)=yb(i);   yv(1)=yb(1);
    else
      sv(1)=spar(i-1); sv(2)=spar(i); sv(3)=spar(i+1)
      xv(3)=xb(i-1);   xv(2)=xb(i);   xv(1)=xb(i+1);
      yv(3)=yb(i-1);   yv(2)=yb(i);   yv(1)=yb(i+1);
    endif
    mA(:,2)=sv(:)
    mA(1,3)=sv(1)*sv(1);
    mA(2,3)=sv(2)*sv(2);
    mA(3,3)=sv(3)*sv(3);
    pc = xv; pA = mA
    CALL DGESVX('Equilibration','No transpose',3,1,pA,3,AF,3,ipiv,   &
           equed,mR,mC,pc,3,xx,3,RCOND,FERR,BERR,WORK,IWORK,INFO)
    if (info .ne. 0) then
      print '(1x,12(e13.6,1x))', rcond
      stop
    endif
    mpoly(2:4,i) = xx(:)
!
    pc = yv; pA = mA
    CALL DGESVX('Equilibration','No transpose',3,1,pA,3,AF,3,ipiv,   &
           equed,mR,mC,pc,3,xx,3,RCOND,FERR,BERR,WORK,IWORK,INFO)
    if (info .ne. 0) then
      print '(1x,12(e13.6,1x))', rcond
      stop
    endif
    mpoly(5:7,i) = xx(:)
    if (i.eq.1) then
      pc(1) = two*mpoly(4,i)*two*cpi+mpoly(3,i)
      pc(2) = two*mpoly(7,i)*two*cpi+mpoly(6,i)
    else
      pc(1) = two*mpoly(4,i)*mpoly(1,i)+mpoly(3,i)
      pc(2) = two*mpoly(7,i)*mpoly(1,i)+mpoly(6,i)
    endif
    tmp=sqrt(pc(1)*pc(1)+pc(2)*pc(2))
    nvx(i)=-pc(2)/tmp; nvy(i) = pc(1)/tmp
  enddo
!
  return
end subroutine getNormal
!
subroutine InitLocCell(xpt, ypt, fsx,fsy, jv)
  integer, parameter :: arc_samples=10000
  double precision, dimension(mpts) :: xpt, ypt, fsx,fsy,jv
!  double precision, dimension(mcoor), intent(out) :: 
!
  double precision :: r, rb, ra, ang, ds, xcmp, ycmp, xcma, ycma, shape_a
  double precision :: arc_theta(0:arc_samples),arc_length(0:arc_samples)
  double precision :: dtheta,target_arc,arc_fraction
  logical :: iReadIn
  integer :: i, j, ient, jent, il, n, sample_index
  double precision :: pa, pb, pc, theta
!
!
  pi = cpi
  ang= two*pi/dble(nring)
  ds = dl/sin(half*ang)
  r  = initial_cell_radius
  ds = cpi*two/dble(nring)
  rb = r/sqrt(initial_cell_axis_ratio)
  ra = r*sqrt(initial_cell_axis_ratio)
  shape_a=initial_cell_shape_factor
  ! Uniform-angle points crowd at the major-axis tips of this ellipse and can
  ! leave multiple markers in one cut cell.  Build a fine cumulative-arclength
  ! table so the requested same-area ellipse retains spacing comparable to hg.
  dtheta=two*pi/real(arc_samples,kind=dp)
  arc_theta(0)=zero
  arc_length(0)=zero
  do i=1,arc_samples
    theta=(real(i,kind=dp)-half)*dtheta
    arc_theta(i)=real(i,kind=dp)*dtheta
    if(abs(shape_a)>epsilon(one))then
      ! Li--Yao--Mori--Sun Fig. 4 material-shape family.  The added
      ! cos(theta)^2 term produces front--rear asymmetry while preserving
      ! area pi*r^2 and transverse width 2*r.
      arc_length(i)=arc_length(i-1)+r*sqrt( &
           (sin(theta)+two*shape_a*cos(theta)*sin(theta))**2+ &
           cos(theta)**2)*dtheta
    else
      arc_length(i)=arc_length(i-1)+sqrt((ra*sin(theta))**2+ &
           (rb*cos(theta))**2)*dtheta
    endif
  enddo
  pa = 9.0d0
  pb = 160.d0
  pc = 3.0d0
!!  pa = 4.0d0
!!  pb = 64.d0
!!  pc = 3.0d0

  ctild=-dl*dl*sin(ang)
  !!ctild=0. ! no preferential curvature
  cs(1) = ctild

  !!xcma=0.35d0
  xcma=half
  ycma=half
  !!ycma=0.25d0
  !!ycma=0.125d0
!
  do n = 1, mpts
    il  = (n-1)/nring + 1 ! index of IB objects
    if (il < 3) then
!    xcmp = xcma + dble(il-1)*2.5d0*ra;
!    ycmp = ycma
    xcmp = half -1.5d0*ra + dble(il-1)*3.1d0*ra;
    ycmp = half
    else
    xcmp = half 
    ycmp = half - 1.50*ra + dble(il-3)*3.1d0*ra
    endif
    ! Marker coordinates are absolute.  On [-1,1]^2 the midpoint is zero;
    ! half*(xmax-xmin) alone would incorrectly place the center at x=1.
    xcmp = xmid+initial_cell_center_offset_x
    ycmp = ymid+initial_cell_center_offset_y
! 
    i   = (n-2*nfil) - (n-2*nfil-1)/nring*nring ! local pointer in each IB object
    target_arc=real(i-1,kind=dp)*arc_length(arc_samples)/ &
         real(nring,kind=dp)
    sample_index=1
    do while(arc_length(sample_index)<target_arc)
      sample_index=sample_index+1
    enddo
    arc_fraction=(target_arc-arc_length(sample_index-1))/ &
         (arc_length(sample_index)-arc_length(sample_index-1))
    theta=arc_theta(sample_index-1)+arc_fraction*dtheta
    if(abs(shape_a)>epsilon(one))then
      xpt(n) = xcmp + r*(cos(theta)+shape_a*cos(theta)**2)
      ypt(n) = ycmp + r*sin(theta)
    else
      xpt(n) = xcmp + ra*cos(theta)
      ypt(n) = ycmp + rb*sin(theta)
    endif
      !!xpt(n) = xcmp - cpi*(pa+cos(pc*theta))*cos(theta)/pb
      !!ypt(n) = ycmp - cpi*(pa+cos(pc*theta))*sin(theta)/pb
  enddo ! end loop over entire list of IB points

!  iReadIn = .true.
!
!  if (iReadIn) then
!    open(unit = 100, file ='ini_cell.txt', status = 'old', action = 'read')
!    do i = 1, mpts
!      read(100,*)xpt(i), ypt(i)
!      print '(1x, i5, 1x, 2(e22.14,1x))', i, xpt(i), ypt(i)
!    enddo
!    close(100)
!
!    !!dbg stop
!  endif

  oxpt = xpt; oypt = ypt; ! (xpt,ypt) are global vars defined in parameters.f90
!  xbp = xpt; ybp = ypt;
!  fs = 0.d0; cfs = 0.d0
  call IBForceVec(xpt,ypt,fsx,fsy)
  call IBJacobian(xpt,ypt, jv)
!  call elasticforce(xpt,ypt,fs,cfs,gp) ! gp stores reciprical of jacobian @ IB pts
  !!dbg fs = zero
  !!dbg cfs = zero
!!prt  print *, 'the force?', maxval(fs), maxval(cfs), maxval(gp)
!
  return
end subroutine InitLocCell
!========================================================================
double precision function ddel(r,isel)
!> Discrete Delta function, 1D
!> @param r is scaled (by 1/h) distance from IB pt
!> 
  double precision :: r
  integer :: isel
!
  select case (isel)
  case (0) ! 0<r<1
    ddel = 0.125D0*(3.d0-2.d0*dabs(r)+dsqrt(1.d0+4.d0*dabs(r)-4.d0*r*r))
  case (1) ! 1<r<2
    ddel = 0.125D0*(5.d0-2.d0*dabs(r)-dsqrt(-7.d0+12.d0*dabs(r)-4.d0*r*r))
  case default
    print *, 'No such choice! stop'
    stop
  end select
!
  return
end function ddel
!========================================================================
subroutine newSpread(xb,yb,fx,fy,ftx,fty)
!> @param (xb,yb) are locs of current IB pts,
!> @param (fx,fy) store forces @ IB pts (xb,yb), in order of list pts on each object
!> @param (ftx,fty) are Eulerian variable for output
  implicit none
  double precision, dimension(mpts), intent(in) :: xb,yb, fx,fy
  double precision, dimension(-1:nx+1,-1:ny+1), intent(out) :: ftx,fty
! 
  integer :: i,j,k, n!, il, nlp, nil, njl, info
  integer :: jmarrLR,kmarrLR
  integer :: jmarrTB,kmarrTB
!  integer :: iflg
  integer :: jmp,kmp
  integer :: ix, jx
!
  double precision :: xxsft, yysft, ds, f1, f2, xper100, yper100
  double precision :: x, y, x0, y0, x1, y1, x2, y2
  double precision :: dxl, dyl, rxl, ryl, dxt, dyt, rxt, ryt
  double precision :: uwtLR(4), vwtLR(4), uwtTB(4), vwtTB(4)
!  double precision :: s0, tmp, tp1, tp2
!
  xamin = xmin
  yamin = ymin
  xper100 = (xmax-xmin)*100.0d0
  yper100 = (ymax-ymin)*100.0d0
!
!========================================================================
  ftx = zero; fty = zero
!
  ds=cpi*two/dble(nring)
  do n=1,mpts
    x=xb(n)+xper100
    y=yb(n)+yper100
    x0=xb(n); y0=yb(n);  ! current IB pt
    !
    f1= ds*(fx(n)) ! force for current IB pt along x
    f2= ds*(fy(n)) ! force for current IB pt along y

! LR edges, for MAC configuration
    xxsft = zero
    yysft = half
    jmarrLR=int((x-xamin)*rh+xxsft)
    kmarrLR=int((y-yamin)*rh+yysft)
!
! compute weights
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x1 = xamin+h*(dble(ix)-xxsft); y1 = yamin+h*(dble(jx)-yysft);
    dxl= x1-x0; dyl = y1-y0; rxl= dxl*rh; ryl = dyl*rh
!
    uwtLR(1)=ddel(rxl-1.d0,1); vwtLR(1)=ddel(ryl-1.d0,1);
    uwtLR(2)=ddel(rxl     ,0); vwtLR(2)=ddel(ryl     ,0);
    uwtLR(3)=ddel(rxl+1.d0,0); vwtLR(3)=ddel(ryl+1.d0,0);
    uwtLR(4)=ddel(rxl+2.d0,1); vwtLR(4)=ddel(ryl+2.d0,1);
    uwtLR = uwtLR*rh;          vwtLR = vwtLR*rh;

! for TB edges
    xxsft = half
    yysft = zero
    jmarrTB=int((x-xamin)*rh+xxsft)
    kmarrTB=int((y-yamin)*rh+yysft)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x2 = xamin+h*(dble(ix)-xxsft); y2 = yamin+h*(dble(jx)-yysft);
    dxt= x2-x0; dyt = y2-y0; rxt= dxt*rh; ryt = dyt*rh
!
    uwtTB(1)=ddel(rxt-1.d0,1); vwtTB(1)=ddel(ryt-1.d0,1);
    uwtTB(2)=ddel(rxt     ,0); vwtTB(2)=ddel(ryt     ,0);
    uwtTB(3)=ddel(rxt+1.d0,0); vwtTB(3)=ddel(ryt+1.d0,0);
    uwtTB(4)=ddel(rxt+2.d0,1); vwtTB(4)=ddel(ryt+2.d0,1);
    uwtTB = uwtTB*rh;          vwtTB = vwtTB*rh;
!
!  forces spread to grid points with indices jmarr-1,..,jmarr+2 and kmarr-1,..,kmarr+2.
!
    do k=1,4
    do j=1,4
      jmp=mod((jmarrLR-2+j)+100*nx,nx)
      kmp=mod((kmarrLR-2+k)+100*ny,ny)
      ftx(jmp,kmp) = ftx(jmp,kmp) + uwtLR(j)*vwtLR(k)*f1
      jmp=mod((jmarrTB-2+j)+100*nx,nx)
      kmp=mod((kmarrTB-2+k)+100*ny,ny)
      fty(jmp,kmp) = fty(jmp,kmp) + uwtTB(j)*vwtTB(k)*f2
    enddo
    enddo
  enddo        ! end loop over entire list of IB pts 
!
! now extend periodically for Eulerian grid
!
  do j=1,ny ! 
    ftx(nx  ,j)=ftx(0,j)
    ftx(nx+1,j)=ftx(1,j)
    fty(nx  ,j)=fty(0,j)
    fty(nx+1,j)=fty(1,j)
  enddo
  fty(nx  ,ny+1)=fty(0,ny+1)
  fty(nx+1,ny+1)=fty(1,ny+1)
!
  return
end subroutine newSpread
!========================================================================
!
!=======================================================================
!=======================================================================
!=======================================================================
!========================================================================
subroutine IBspreadS(xpt,ypt,gp,fs,ff)
!> @param fs store forces @ IB pts, in order (f1_x,f1_y,f2_x,f2_y,...)
!> @param gp store Jacobian @ each IB pt
!> 
  implicit none
  double precision, dimension(mcoor), intent(in) :: fs
  double precision, dimension(mpts), intent(in) :: xpt, ypt, gp
  double precision, intent(out) :: ff(-1:nx+1,-1:ny+1,2) ! force on Eulerian grid
  !
!  double precision :: f(-1:nxp2,-1:nyp2,2)
  double precision :: g(-1:nx+1,-1:ny+1,2)
!
! for mac grid 
!
! The spreading is done to an extended domain with 2 extra cells above
! and below the tube walls.  This force will be mapped (with reflection)
! to a smaller force array before being passed to fluidstep.
!
! spreads to 4x4 portion of grid directly
!
! lower left of the spreading grid  is (xamin,yamin)
! upper right of the spreading grid is (xamax,yamax)
!
! ib point should be located between
!    xamin + 2h and xamax-2h
!    yamin + 2h and yamax-2h
! 
  integer :: i,j,k, n, il, nlp, nil, njl, info
  integer :: jmarrLR,kmarrLR
  integer :: jmarrTB,kmarrTB
  integer :: iflg
  integer :: jmp,kmp
  integer :: ix, jx
!
  double precision :: xxsft, yysft, ds, f1, f2, xper100, yper100
  double precision :: x, y, x0, y0, x1, y1, x2, y2
  double precision :: dxl, dyl, rxl, ryl, dxt, dyt, rxt, ryt
  double precision :: uwtLR(4), vwtLR(4), uwtTB(4), vwtTB(4)
  double precision :: s0, tmp, tp1, tp2
!
  xamin = xmin
  yamin = ymin
  xper100 = (xmax-xmin)*100.0d0
  yper100 = (ymax-ymin)*100.0d0
!
!========================================================================
  g(0:nx,0:ny,1:2)=zero
!
  ds=cpi*two/dble(nring)
  do n=1,mpts
    x=xpt(n)+xper100
    y=ypt(n)+yper100
    x0=xpt(n); y0=ypt(n);  ! current IB pt
    !
    il  = (n-2*nfil-1)/nring + 1 ! # of IB object for current IB pt
    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring ! local ptr in current IB object
    nil = nlp + (il-1)*2*nring ! ptr of storage for the x-component of current IB
    njl = nil + nring ! ptr of storage for the y-component of current IB
!
    f1= ds*(fs(2*n-1)) ! force for current IB pt along x
    f2= ds*(fs(2*n  )) ! force for current IB pt along y

! LR edges, for MAC configuration
    xxsft = zero
    yysft = half
    jmarrLR=int((x-xamin)*rh+xxsft)
    kmarrLR=int((y-yamin)*rh+yysft)
!
! compute weights
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x1 = xamin+h*(dble(ix)-xxsft); y1 = yamin+h*(dble(jx)-yysft);
    dxl= x1-x0; dyl = y1-y0; rxl= dxl*rh; ryl = dyl*rh
!
    uwtLR(1)=ddel(rxl-1.d0,1); vwtLR(1)=ddel(ryl-1.d0,1);
    uwtLR(2)=ddel(rxl     ,0); vwtLR(2)=ddel(ryl     ,0);
    uwtLR(3)=ddel(rxl+1.d0,0); vwtLR(3)=ddel(ryl+1.d0,0);
    uwtLR(4)=ddel(rxl+2.d0,1); vwtLR(4)=ddel(ryl+2.d0,1);
    uwtLR = uwtLR*rh;          vwtLR = vwtLR*rh;

! for TB edges
    xxsft = half
    yysft = zero
    jmarrTB=int((x-xamin)*rh+xxsft)
    kmarrTB=int((y-yamin)*rh+yysft)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x2 = xamin+h*(dble(ix)-xxsft); y2 = yamin+h*(dble(jx)-yysft);
    dxt= x2-x0; dyt = y2-y0; rxt= dxt*rh; ryt = dyt*rh
!
    uwtTB(1)=ddel(rxt-1.d0,1); vwtTB(1)=ddel(ryt-1.d0,1);
    uwtTB(2)=ddel(rxt     ,0); vwtTB(2)=ddel(ryt     ,0);
    uwtTB(3)=ddel(rxt+1.d0,0); vwtTB(3)=ddel(ryt+1.d0,0);
    uwtTB(4)=ddel(rxt+2.d0,1); vwtTB(4)=ddel(ryt+2.d0,1);
    uwtTB = uwtTB*rh;          vwtTB = vwtTB*rh;
!
!!deb      print '(1x,"which part ", 2(i5,1x),1024(e14.6,1x))', jmarrTB, kmarrTB, &
!!deb        x,y,xpt(n),ypt(n), xper100, yper100, xxsft,yysft, uwtLR(2), vwtLR(2), f1, f2
!  spread these forces to Eulerian grids (MAC configuration) now
!
!  forces spread to grid points with indices jmarr-1,..,jmarr+2 and kmarr-1,..,kmarr+2.
!
    do k=1,4
    do j=1,4
      jmp=mod((jmarrLR-2+j)+100*nx,nx)
      kmp=mod((kmarrLR-2+k)+100*ny,ny)
      g(jmp,kmp,1) = g(jmp,kmp,1) + uwtLR(j)*vwtLR(k)*f1
      jmp=mod((jmarrTB-2+j)+100*nx,nx)
      kmp=mod((kmarrTB-2+k)+100*ny,ny)
      g(jmp,kmp,2) = g(jmp,kmp,2) + uwtTB(j)*vwtTB(k)*f2
    enddo
    enddo
  enddo        ! end loop over entire list of IB pts 
!
!!  print '(1x,"Difference in force: ", 6(e16.8,1x))', &
!!   & norm2(f(0:nx,0:ny,1)-g(0:nx,0:ny,1)), norm2(f(0:nx,0:ny,2)-g(0:nx,0:ny,2)), &
!!   & norm2(f(0:nx,0:ny,1)), norm2(g(0:nx,0:ny,1)),  &
!!   & norm2(f(0:nx,0:ny,2)), norm2(g(0:nx,0:ny,2)) 
!
! now extend periodically for Eulerian grid
!
  do j=1,ny ! 
    g(nx  ,j,1)=g(0,j,1)
    g(nx+1,j,1)=g(1,j,1)
    g(nx  ,j,2)=g(0,j,2)
    g(nx+1,j,2)=g(1,j,2)
  enddo
  g(nx  ,ny+1,2)=g(0,ny+1,2)
  g(nx+1,ny+1,2)=g(1,ny+1,2)
!
  ff(-1:nx+1,-1:ny+1,1:2) = g(-1:nx+1,-1:ny+1,1:2)
!========================================================================
! Delta function is new
!========================================================================
!
  return
end subroutine IBspreadS
!========================================================================
!========================================================================
subroutine newInterpS(xb,yb,u0,v0,tx,ty)
!  interpolate Eulerian variable (u0,v0) defined on MAC grids to the
!  Lagrange variable tgt, which is defined @ IB locations (xb,yb). 
  implicit none
  double precision, dimension(mpts), intent(in):: xb, yb
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in) :: u0, v0
  double precision, dimension(mpts), intent(out) :: tx,ty !<-this is the output
!
  integer :: jmarrLR,kmarrLR
  integer :: jmarrTB,kmarrTB
  integer :: n,ngrid, nlp, nil, njl, il, info, iflag
  integer :: j,k, iflg
  integer :: jmp,kmp
  double precision :: xper100,yper100
  double precision :: uwtLR(4),vwtLR(4)
  double precision :: uwtTB(4),vwtTB(4)
  double precision :: x,y, x0,y0, x2,y2, x1, y1, xp, yp, xxsft, yysft
  double precision :: dxl, dyl, rxl, ryl, dxt, dyt, rxt, ryt
  integer :: ix, jx
!
  xamin = xmin; yamin = ymin
  xper100 = (xmax-xmin)*100.0d0
  yper100 = (ymax-ymin)*100.0d0
!
  do n=2*nfil+1,mpts
    x0 = xb(n); y0 = yb(n)
    x  = x0 + xper100; y = y0 + yper100
    !
! for LR edges
    xxsft = 0.0d0
    yysft = 0.5d0
    jmarrLR=int((x-xamin)*rh+xxsft)
    kmarrLR=int((y-yamin)*rh+yysft)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x1 = xamin+h*(dble(ix)-xxsft); y1 = yamin+h*(dble(jx)-yysft);
    dxl= x1-x0; dyl = y1-y0; rxl= dxl*rh; ryl = dyl*rh
!
    uwtLR(1)=ddel(rxl-1.d0,1); vwtLR(1)=ddel(ryl-1.d0,1);
    uwtLR(2)=ddel(rxl     ,0); vwtLR(2)=ddel(ryl     ,0);
    uwtLR(3)=ddel(rxl+1.d0,0); vwtLR(3)=ddel(ryl+1.d0,0);
    uwtLR(4)=ddel(rxl+2.d0,1); vwtLR(4)=ddel(ryl+2.d0,1);
    uwtLR = uwtLR*rh;          vwtLR = vwtLR*rh;

! for TB edges
    xxsft = 0.5d0
    yysft = 0.0d0
    jmarrTB=int((x-xamin)*rh+xxsft)
    kmarrTB=int((y-yamin)*rh+yysft)
!
    ix = int((x0-xamin)*rh+xxsft); jx = int((y0-yamin)*rh+yysft)
    x2 = xamin+h*(dble(ix)-xxsft); y2 = yamin+h*(dble(jx)-yysft);
    dxt= x2-x0; dyt = y2-y0; rxt= dxt*rh; ryt = dyt*rh
!
    uwtTB(1)=ddel(rxt-1.d0,1); vwtTB(1)=ddel(ryt-1.d0,1);
    uwtTB(2)=ddel(rxt     ,0); vwtTB(2)=ddel(ryt     ,0);
    uwtTB(3)=ddel(rxt+1.d0,0); vwtTB(3)=ddel(ryt+1.d0,0);
    uwtTB(4)=ddel(rxt+2.d0,1); vwtTB(4)=ddel(ryt+2.d0,1);
    uwtTB = uwtTB*rh;          vwtTB = vwtTB*rh;
!
!compute weights
!
!  interpolate these Lagrangian quantity by computing
!  the weighted sum that gives the values at (xb,yb).
!
    tx(n)=0.d0
    ty(n)=0.d0
!
    do k=1,4
    do j=1,4
      jmp=mod((jmarrLR-2+j)+100*nx,nx)
      kmp=mod((kmarrLR-2+k)+100*ny,ny)
      !!CosDelt tx(n)=tx(n)+xwtLR(j)*ywtLR(k)*u0(jmp,kmp) 
      tx(n)=tx(n)+uwtLR(j)*vwtLR(k)*u0(jmp,kmp) 
      jmp=mod((jmarrTB-2+j)+100*nx,nx)
      kmp=mod((kmarrTB-2+k)+100*ny,ny)
      !!CosDelt ty(n)=ty(n)+xwtTB(j)*ywtTB(k)*v0(jmp,kmp) 
      ty(n)=ty(n)+uwtTB(j)*vwtTB(k)*v0(jmp,kmp) 
    enddo
    enddo
!    tgt(n,:)=tgt(n,:)*h*h
  enddo
!  tx = tx*hg*hg;
!  ty = ty*hg*hg
!
  return
end subroutine newInterpS

!========================================================================
!========================================================================
!========================================================================
!
! NOTE: petmove (a backward-Euler, single-solve implicit dX/dt update) used
! to live here. It was dead code (only reachable via a commented-out call
! in fmain.f90) and has been retired in favor of the actual running
! algorithm, now in fsisolve.f90 (AdvanceFSI/zMatMul: midpoint-rule,
! Picard-iterated, full fluid+elastic coupling).
!=======================================================================
subroutine normalproject(fx,fy, nvx,nvy, px,py)
!> @param (nvx,nvy) normal vector
!> @param (fx,fy) vector field 
!> @param (ps,py) output
  double precision, dimension(mpts), intent(out) :: px,py
  double precision, dimension(mpts), intent(in)  :: fx,fy, nvx,nvy
!
  integer :: j
  double precision :: tmp

  do j = 1, mpts
    tmp = fx(j)*nvx(j)+fy(j)*nvy(j)
    px(j) = tmp*nvx(j)
    py(j) = tmp*nvy(j)
  enddo
!
  return
end subroutine normalproject
!
subroutine hadamard(vx,vy, vo)
  double precision, dimension(mpts), intent(out) :: vo
  double precision, dimension(mpts), intent(in)  :: vx, vy
!
  integer :: j
!
  do j = 1, mpts
    vo(j) = vx(j)*vy(j)
  enddo
!
  return
end subroutine hadamard
!
!=======================================================================
! NOTE: linMatMul (petmove's shell-matrix operator) used to live here;
! retired alongside petmove -- see fsisolve.f90's zMatMul instead.
!=======================================================================
double precision function evalmpoly(s, n, mk, idir, icase, info)
! evaluate function values along IB points, by using interpolation stored
! in mk(1:7,1:nring), 
! idir: 0<-x; 1<-y
! icase: 0<- original interpolation; 1<- 1st derivative; 2<- 2nd derivative
  integer :: lda, n
  double precision :: s, mk(7,n)
  integer :: idir, icase, info
!
  double precision :: ck, tupi, eps
  integer :: i, k, inc, ijump
! 
  tupi = 2.d0*cpi
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
    if (ck < mk(1,k)) then
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
      if ( ck .ge. mk(1,k) .and. ck < mk(1,k+1)) then
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
      if ( ck .ge. mk(1,k) .and. ck < tupi) then
        ijump = 1
      else
        print *, "double precisionly?", ck, s, mk(1,1)
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
    evalmpoly = mk(inc,k) + (mk(inc+1,k) + mk(inc+2,k)*ck)*ck
  case (1)
!!  if (icase .eq. 1) then
    !!evalmpoly = mk(inc,k,il) + mk(inc+1,k,il)*ck + mk(inc+2,k,il)*ck*ck
    evalmpoly =  mk(inc+1,k) + 2.d0*mk(inc+2,k)*ck
!!  else if (icase .eq. -1) then ! linear interpolation
!!  endif
  !!prt print *, 'ck = ', ck, ' k = ', k, 'inc = ', inc, 'idir = ', idir, 's= ', s
  case (2)
    evalmpoly = 2.d0*mk(inc+2,k)
  case default
    evalmpoly = mk(inc,k) + (mk(inc+1,k) + mk(inc+2,k)*ck)*ck
    print *, 'no such option, evaluate original interpolation'
  end select
  if (ijump .eq. 1) then
    info = 0
  else
    print '(1x,"why here?", 4(i3,1x), 1024(e16.8,1x))',k, info, n, ijump, s, ck, evalmpoly
    stop
  endif
!
  return
end function evalmpoly
!
subroutine ResetPins(xb,yb, pinu,pinv, mk, nb)
! ResetPins: reset grid pin for in/out cell interface which (xb,yb) moves 
!   so Brinkman's term can be computed correctly (discontinous constants)
! on input, id & idf store info from previous time step n (before sub iteration step)
! on input, mk & nb store bdry interpolation and normal respectively 
! on input, isel = -1 for initial setup, when oid & kdf do not contain info(NOT USIng)
! on output, pinu & pinv store current tagging info, which may be the same
!   as (jdu,jdv)
!-----B--1----+----2----+----3----+----4----+----5----+----6----+----7-E
  implicit none
!  integer, dimension(-1:nx+1,-1:ny+1) :: id, idf
  integer, dimension(-1:nx+1,-1:ny+1) :: pinu, pinv
  integer :: isel
  double precision :: mk(7,nring), nb(nring,2)
  double precision, dimension(mpts), intent(in) :: xb,yb
!  type(ibpt), dimension(:) :: ibary(nent-2)
!
  logical :: is_intra, all_in(4)
  integer :: i, j, k, is, js, kk, info, il, il_save, ijump, jjump, ip, jp
  double precision :: x, y, dx, dy, xt, yt, vx, vy, x0, y0
  double precision :: tmp, tp1, tp2, tp3, tp4, s0, s1, tpv(nib)
  double precision :: ushift, vshift
!
!  print*, 'are we here in ResetPins: mk1'
  dx = hg; dy = hg
  pinu = jdu; pinv = jdv ! (jdu,jdv) global vars for Pin

!  idu = 0; idv = 0; jdu = 0; jdv = 0
  ushift = zero; vshift = half
  do j = 2, ny-1
  do i = 1, nx-1
    if (dble(idf(i,j)+idf(i+1,j)) > -0.5) then
      if (dble(idf(i,j) + idf(i+1,j)) >1.1 .and. dble(idf(i,j+1)+idf(i+1,j+1))>1.1 .and. dble(idf(i,j-1)+idf(i+1,j-1))>1.1) then
        pinu(i,j) = 1 !jdu
        exit
      endif
      x = xmin+dx*(dble(i) - ushift); y = ymin+dy*(dble(j)-vshift)

      do il = 1, nib
        tmp = 10.d0; tp3 = 10.d0; kk = 1
!
        do k = 1, nring
          !tp1 = ((x-curr%x)**2.d0+(y-curr%y)**2.d0)
          tp1 = ((x-xb(k))**2.d0+(y-yb(k))**2.d0)

!          print  '(1x,"dist", e14.6," tmp",e14.6, " kk",i3, " x,y",       &
!                             2(e14.6,1x))',tp1, tmp, kk
          if (tp1 <= tmp) then
            kk = k
            tmp = tp1
            !xt = curr%x; yt = curr%y; s0 = mk(1,kk,il); s1 = s0
            xt = xb(k); yt = yb(k); s0 = mk(1,kk); s1 = s0
          endif
!          curr => curr%prev
        enddo

        tmp = (x-xt)*nb(kk,1)+(y-yt)*nb(kk,2)

        if (tmp < 0. ) then
!          id(i,j) = 1 ! inside of platelet
          all_in(il) = .true.      
        else
!          id(i,j) = 0
          all_in(il) = .false.
        endif
      enddo ! end of sweep on platelets
!
      is_intra = .false. 
      ijump = 0
      do il = 1, nib
        is_intra = is_intra .or. all_in(il)
      enddo
      if (is_intra) then
        pinu(i,j) = 1
      endif
    endif
  enddo
  enddo
  ushift = half; vshift = zero
  do j = 2, ny-1
  do i = 2, nx-1
    if (idf(i,j) + idf(i,j+1) >-0.1) then
      if (dble(idf(i,j) + idf(i,j+1)) >1.1 .and. dble(idf(i+1,j)+idf(i+1,j+1))>1.1 .and. dble(idf(i-1,j)+idf(i-1,j+1))>1.1) then
        pinv(i,j) = 1
        exit
      endif
      x = xmin+dx*(dble(i) - ushift); y = ymin+dy*(dble(j)-vshift)

      do il = 1, nib
        tmp = 10.d0; tp3 = 10.d0; kk = 1
!        curr => ibary(il)%p
!
        do k = 1, nring
          !tp1 = ((x-curr%x)**2.d0+(y-curr%y)**2.d0)
          tp1 = ((x-xb(k))**2.d0+(y-yb(k))**2.d0)
          if (tp1 <= tmp) then
            kk = k
            tmp = tp1
            !xt = curr%x; yt = curr%y; s0 = mk(1,kk,il); s1 = s0
            xt = xb(k); yt = yb(k); s0 = mk(1,kk); s1 = s0
          endif
!          curr => curr%prev
        enddo
        tmp = (x-xt)*nb(kk,1)+(y-yt)*nb(kk,2)
        if (tmp < 0. ) then
!          id(i,j) = 1 ! inside of platelet
          all_in(il) = .true.      
        else
!          id(i,j) = 0
          all_in(il) = .false.
        endif
      enddo ! end of sweep on platelets
!
      is_intra = .false. 
      do il = 1, nib
        is_intra = is_intra .or. all_in(il)
      enddo
      if (is_intra) then
        pinv(i,j) = 1
      endif
    endif
  enddo
  enddo
!========================================================================
! End of identifying fluid/solid points
!========================================================================
!  print *, 'Are we done in tagging: mk2'

  return
end subroutine ResetPins
!========================================================================
double precision function pusher(s)
  implicit none
! 
  double precision :: s
!
  double precision :: shs, shw, shl, sts, stw, stl, sul, suw, sus, sbl, sbw, sbs
!p
  shs=-.0080d0;   sts=-1.0*shs
  sus=-shs*0.0;  sbs= 1.0*sus
  shl=cpi*.0 ;   stl=cpi*1.0
  sul=cpi*.80;   sbl=cpi*1.20
  shw=cpi*.21;   stw=shw*1.0;
  suw=shw*.5 ;   sbw=suw*1.

!!o  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+0. *exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!!o      sts*exp(-.5*(s-stl)**2./stw**2.)
!working  tmp=shs*exp(-.5*(s-shl)**2./shw**2.)+shs*exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!working      sts*exp(-.5*(s-stl)**2./stw**2.)
!  pusher= shs*exp(-.5*(s-shl)**2./shw**2.)+shs*exp(-.5*(s-shl-tupi)**2./shw**2.)+&
!      sts*exp(-.5*(s-stl)**2./stw**2.)+ &
!      sus*exp(-.5*(s-sul)**2./suw**2.)+sbs*exp(-.5*(s-sbl)**2./sbw**2.)
  pusher=shs*tanh(half*cos(s-shl))
!  pusher = 0.001d0
  return
end function pusher
!========================================================================
subroutine curvePoly(mk,conv, convsq)
  double precision, dimension(:,:,:) :: mk(7,nring)
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

!  do il = 1, nib
    do i = 1, nring-1
      s0  = 0.5d0*(mk(1,i)+mk(1,i+1)); 
      dst = (mk(1,i+1)-mk(1,i))
      xp  = evalmpoly(s0, nring, mk, 0, 1, info)
      yp  = evalmpoly(s0, nring, mk, 1, 1, info)
      xpp = evalmpoly(s0, nring, mk, 0, 2, info)
      ypp = evalmpoly(s0, nring, mk, 1, 2, info)
      t1  = (xp*xp+yp*yp)
      cv  = (xp*ypp - yp*xpp)/t1*dst
      tcv = tcv + cv
      conv= conv + abs(cv)
      convsq = convsq + cv*cv
    enddo
    i = nring
    s0 = 0.5d0*(mk(1,i)+tupi); ! note, we have s0< tupi here
    dst = (tupi-mk(1,i))
    xp  = evalmpoly(s0, nring, mk, 0, 1, info)
    yp  = evalmpoly(s0, nring, mk, 1, 1, info)
    xpp = evalmpoly(s0, nring, mk, 0, 2, info)
    ypp = evalmpoly(s0, nring, mk, 1, 2, info)
    t1  = (xp*xp+yp*yp)
    cv  = (xp*ypp - yp*xpp)/t1*dst
    tcv = tcv + cv
    conv= conv + abs(cv)
    convsq = convsq + cv*cv
!  enddo
  
!!  print '(1x,"absolute CURVATURE: ", e12.6, "curvature: ", 10(e12.6,1x))', conv/tupi, tcv/tupi, mk(1,1,il), mk(1,nring,il)
!
  return
end subroutine curvePoly
!
end module IBmod
