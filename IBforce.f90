!> Computing forces @ each IB pts from IB pt locations (xpt,ypt)
!> 
!> @author Lingxing Yao
!> @version 0.1
!> 
module IBforce
  use, intrinsic :: iso_c_binding
  use parameters

implicit none

  double precision, parameter :: h=hg

  private

  logical,parameter :: consttension = .false.
  double precision, parameter :: dls = two*cpi/dble(nring), rds=one/dls

  !public :: elasticforce, IBJacobXY, IBJacobian, IBForceVec
  public ::  IBJacobXY, IBJacobian, IBForceVec
  ! IBJacobian return |dX/dsigma| @ IB pts
  ! IBJacobXY return J(X)*Y, with Jacobian of F(X) being J(X)=(dF(X)/dX)
  ! IBForceVec return force vectors: two vectors for x & y direction
  ! elasticforce return force vector in single vector form and 2 components of force @ IB pts stored continuously

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!tmpsubroutine elasticforce(xb,yb,fs,cfs,gp)
!tmp!> Computing force at each IB pt
!tmp!> 
!tmp!> @param xb  x-coor @ IB pts, mpts = cent*nring
!tmp!> @param yb  y-coor @ IB pts, mpts = cent*nring
!tmp!> @param fs  calculated force vector @ each IB pts, fs(mcoor = mpts*2), and forces at each IB pt are 
!tmp!>   stored continuously with (f1_x,f1_y, f2_x,f2_y,.....)
!tmp!> @param cfs  scaled force vector @ each IB pts, mcoor = mpts*2, cfs=fs*gp, same storage arrangement as in fs
!tmp!> @param gp  jacobian at each IB pt, in shape (mpts)
!tmp
!tmp!integer, intent(in) :: isel
!tmpinteger      :: j
!tmpdouble precision  :: xb(mpts) ! x-coor @ IB pts, mpts = cent*nring
!tmpdouble precision  :: yb(mpts)
!tmpdouble precision  :: fs(mcoor) ! force vector @ IB pts, mcoor = mpts*2
!tmpdouble precision  :: cfs(mcoor) ! scaled force vector=force*Jacobian = fs*gp 
!tmpdouble precision  :: g(mpts,2) ! force vector in (fx,fy) form
!tmpdouble precision  :: gp(mpts) ! Jacobian 
!tmp!
!tmp! now we evaluate force density at IB locations:
!tmp!
!tmp!   write(*,*)' elasticforce'
!tmp      
!tmp      call ngrad(xb,yb,g,gp) ! gp stores Jacobian @ each IB pt
!tmp!      write(*,*)' after grad max g=',maxval(abs(g))
!tmp!
!tmp! calculate entity force density (force/length)
!tmp!
!tmp      cfs = 0.
!tmp      do j=1,mpts
!tmp         fs(2*j-1)=-g(j,1)
!tmp         fs(2*j)  =-g(j,2)
!tmp      enddo
!tmp      !!oldIdx do j=2*nfil+1,npts
!tmp      do j=1,mpts ! compute f*|dx/ds|
!tmp         cfs(2*j-1)=-g(j,1)*gp(j)
!tmp         cfs(2*j)  =-g(j,2)*gp(j)
!tmp      enddo
!tmp!      write(*,*)' IB obj max fs=',maxval(abs(fs))
!tmp!!prt      print *, 'force max ', maxval(gp), isel
!tmp  return
!tmpend subroutine elasticforce
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!123456789212345678931234567894123456789512345678961234567897123456789
subroutine ngrad(x,y,g,gp)
!> computing force
!> @author Lingxing Yao
!> @param x input: x-coord @ each IB pt, x(mpts)
!> @param y input: y-coord @ each IB pt, y(mpts)
!> @param g output: calculated force, g(mpts,2)
!> @param gp output: Jacobian factor @ each IB pt, gp(mpts)
!
  double precision, dimension(mpts), intent(in) :: x,y
  double precision, dimension(mpts), intent(out) :: gp
  double precision, intent(out)  :: g(mpts,2)
  !
  integer :: i,im,ip, n, il, nlp!, nil, njl
!  double precision  :: ecv,g1,g2,gp1,gp2
  double precision  :: sm,xm1,xp,xp1, sp
  double precision  :: ym1,yp,yp1, gd
!  double precision  :: xs1,xq,xa1, ys1,yq,ya1
!
  double precision  :: tmp
  double precision, dimension(mpts) :: xx, yy, cx, cy
!
  g = zero
  gp= zero
!
!=======================================================================
!!  find difference of current IB pts and preferred shape
!!  do i = 1, nring
!!    print '(1x,10(e22.16,1x))', x(i), y(i), xbp(i), ybp(i), x(i)-xbp(i), y(i)-ybp(i)
!!  enddo
  
!!prt  print '(1x,10(e22.16,1x))', x(1), y(1), xbp(1), ybp(1), x(1)-xbp(1), y(1)-ybp(1)
!!prt  print *, ''
!!  if (abs(sum(x-xbp))+abs(sum(y-ybp))>1d-10 ) then
!    xx= x-xbp; yy= y-ybp ! using difference between input IB location and preferred shape,
                         ! could be set to xx=x, yy=y so no preferred shape used
                         ! xbp and ybp are initial IB pts location
    xx =x; yy = y; 
!!  else
!!    xx= x; yy= y
!!  endif
!=======================================================================
!
!=======================================================================

!=======================================================================
do n = 1, mpts  ! loop over entire list of IB pts
  il  = (n-2*nfil-1)/nring + 1 ! index of IB objects
  nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring ! local pointer in each IB object
!  nil = nlp + (il-1)*2*nring ! pointer for x-coord of the current (n-th) IB pt
!  njl = nil + nring ! pointer for y-coord of the current (n-th) IB pt
!
  i = nlp
  if (i>1 .and. i< nring) then
    im = n-1; ip = n+1
  elseif (i.eq.1) then
    im = n+nring-1; ip = n+1
  elseif (i.eq.nring)then
    im = n-1; ip = n-nring + 1
  endif
  cx(n)=xx(ip)-two*xx(n)+xx(im)
  cy(n)=yy(ip)-two*yy(n)+yy(im)
enddo !end loop over entire list of  IB pts
cx=cx*rds*rds; cy=cy*rds*rds

tmp = rsl*dls ! rest lenth = rsl

do n = 1, mpts  ! loop over entire list of IB pts
  il  = (n-2*nfil-1)/nring + 1 ! index of IB objects
  nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring ! local pointer in each IB object
!  nil = nlp + (il-1)*2*nring ! pointer for x-coord of the current (n-th) IB pt
!  njl = nil + nring ! pointer for y-coord of the current (n-th) IB pt

  i = nlp
    if (i>1 .and. i< nring) then
      im = n-1; ip = n+1
    elseif (i.eq.1) then
      im = n+nring-1; ip = n+1
    elseif (i.eq.nring)then
      im = n-1; ip = n-nring+1
    endif
    xp=xx(n); xp1=xx(ip); xm1=xx(im)
    yp=yy(n); yp1=yy(ip); ym1=yy(im)
!    xq=x(n); xa1=x(ip); xs1=x(im)
!    yq=y(n); ya1=y(ip); ys1=y(im)
!
    sp=sqrt((xp1-xp)*(xp1-xp) + (yp1-yp)*(yp1-yp))
    sm=sqrt((xm1-xp)*(xm1-xp) + (ym1-yp)*(ym1-yp))
!
!=======================================================================
! w/o rest length
    !g(n,1)=-rds*((one-dls/sp)*rds*(xp1-xp)-(one-dls/sm)*rds*(xp-xm1))
    !g(n,2)=-rds*((one-dls/sp)*rds*(yp1-yp)-(one-dls/sm)*rds*(yp-ym1))
    !g(n,:)=sw(1)*g(n,:)
    !constant tension
    !gd=(sp+sm)/(dls*two)
    !g(n,1)=-sw(il)*(cx(n)/gd-half*(xp1-xm1)*rds/(gd*gd)*(sp-sm)*rds)
    !g(n,2)=-sw(il)*(cy(n)/gd-half*(yp1-ym1)*rds/(gd*gd)*(sp-sm)*rds)
    ! version 2
    g(n,1)=-sw(il)*((xp1-xp)/sp-(xp-xm1)/sm)*rds
    g(n,2)=-sw(il)*((yp1-yp)/sp-(yp-ym1)/sm)*rds
    !!
!=======================================================================
! with bending force
!    g(n,1)=g(n,1)+sb(1)*rds*rds*(cx(im)-two*cx(n)+cx(ip))
!    g(n,2)=g(n,2)+sb(1)*rds*rds*(cy(im)-two*cy(n)+cy(ip))
!=======================================================================
    gp(n)=dls*two/((sm)+(sp))
!=======================================================================
enddo ! end loop over entire IB list
!
  return
end subroutine ngrad ! old elastic force implementation
!=======================================================================
!123456789212345678931234567894123456789512345678961234567897123456789
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================
subroutine IBJacobian(xb,yb, jv)
  double precision, dimension(mpts), intent(in)  :: xb, yb
  double precision, dimension(mpts), intent(out) :: jv
!
  integer :: n, il, nlp, im, ip
!
  double precision :: sm,sp, xm1,ym1, xq,yq, xp1,yp1
!
  il = 1
  do n = 1, mpts  ! loop over entire list of IB pts of all IB objs
!    il  = (n-2*nfil-1)/nring + 1 ! index of IB objects
!    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring ! local pointer of each IB object
    nlp = n

    if (nlp>1 .and. nlp< nring) then
      im = n-1; ip = n+1
    elseif (nlp.eq.1) then
      im = n+nring-1; ip = n+1
    elseif (nlp.eq.nring)then
      im = n-1; ip = n-nring+1
    endif
    xq=xb(n); xp1=xb(ip); xm1=xb(im)
    yq=yb(n); yp1=yb(ip); ym1=yb(im)
  !
    sp=sqrt((xp1-xq)*(xp1-xq) + (yp1-yq)*(yp1-yq))
    sm=sqrt((xm1-xq)*(xm1-xq) + (ym1-yq)*(ym1-yq))
    jv(n) = (sm+sp)/(two*dls)
  enddo
!
  return
end subroutine IBJacobian
!
!=======================================================================
!123456789212345678931234567894123456789512345678961234567897123456789
subroutine IBJacobXY(xb,yb,xd,yd, jvx,jvy)
  !apply Jacob matrix 
  !@(xb,yb) Jacobian base vector; along (xd,yd) direction 
  !the result is in (jvx,jvy)
  double precision, dimension(mpts), intent(in)  :: xb, yb, xd, yd
  double precision, dimension(mpts), intent(out) :: jvx, jvy
!
  double precision, dimension(mpts) ::  jpx, jpy
!
  integer :: il, nlp, ip, im, n
  double precision :: xbq,ybq, xbp,ybp, xbm,ybm, xdq,ydq, xdp,ydp, xdm,ydm, gsc
  double precision :: sp, sm
!
  il = 1
  if (consttension) then! T=k; G=-k
    gsc=-sw(il)
  else ! T(s)=k*(s-s_ref); G(s)=k
    gsc= sw(il)
  endif
  do n = 1, mpts
!    il  = (n-2*nfil-1)/nring + 1 ! index of IB objects
!    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring ! local pointer of each IB object

    nlp = n
    if (nlp>1 .and. nlp< nring) then
      im = n-1; ip = n+1
    elseif (nlp.eq.1) then
      im = n+nring-1; ip = n+1
    elseif (nlp.eq.nring)then
      im = n-1; ip = n-nring+1
    endif
    !??q: current n; ??p: n+1; ??m: n-1
    xbq=xb(n); xbp=xb(ip); xbm=xb(im); ybq=yb(n); ybp=yb(ip); ybm=yb(im)
    xdq=xd(n); xdp=xd(ip); xdm=xd(im); ydq=yd(n); ydp=yd(ip); ydm=yd(im)

    sp=sqrt((xbp-xbq)*(xbp-xbq) + (ybp-ybq)*(ybp-ybq))! sq/ds=|D^-X_{i+1}|
    sm=sqrt((xbm-xbq)*(xbm-xbq) + (ybm-ybq)*(ybm-ybq))! sm/ds=|D^-X_{i}|

    !J^0 part
    if (consttension)then ! T=k
      jvx(n)= sw(il)*((xdp-xdq)/sp-(xdq-xdm)/sm)
      jvy(n)= sw(il)*((ydp-ydq)/sp-(ydq-ydm)/sm)
    else ! T=k*(s-s_ref)
      jvx(n)=sw(il)*((sp/dls-membrane_reference_metric)*(xdp-xdq)/sp- &
           (sm/dls-membrane_reference_metric)*(xdq-xdm)/sm)
      jvy(n)=sw(il)*((sp/dls-membrane_reference_metric)*(ydp-ydq)/sp- &
           (sm/dls-membrane_reference_metric)*(ydq-ydm)/sm)
    endif

    !J^1 part
    jvx(n)=jvx(n)+gsc*( ((xbp-xbq)*(xdp-xdq)+(ybp-ybq)*(ydp-ydq))*(xbp-xbq)/ &
      & (sp*sp*sp) -((xbq-xbm)*(xdq-xdm)+(ybq-ybm)*(ydq-ydm))*(xbq-xbm)/(sm*sm*sm))
    jvy(n)=jvy(n)+gsc*( ((xbp-xbq)*(xdp-xdq)+(ybp-ybq)*(ydp-ydq))*(ybp-ybq)/ &
      & (sp*sp*sp) -((xbq-xbm)*(xdq-xdm)+(ybq-ybm)*(ydq-ydm))*(ybq-ybm)/(sm*sm*sm))

    jvx(n)=jvx(n)*rds
    jvy(n)=jvy(n)*rds
  enddo
!
  return
end subroutine IBJacobXY
!
!=======================================================================
!123456789212345678931234567894123456789512345678961234567897123456789
subroutine IBForceVec(xb,yb, Fvx,Fvy)
  double precision, dimension(mpts), intent(in)  :: xb, yb
  double precision, dimension(mpts), intent(out) :: Fvx, Fvy
!
  integer :: n, nlp, il, im, ip
  double precision :: xp1,yp1, xm1,ym1, xq,yq, sp, sm, cx, cy, gd
!  double precision, dimension(mpts) :: jv
!
!
  il = 1
  do n = 1, mpts  ! loop over entire list of IB pts of all IB objs
!    il  = (n-2*nfil-1)/nring + 1 ! index of IB objects
!    nlp = (n-2*nfil) - (n-2*nfil-1)/nring*nring ! local pointer of each IB object

    nlp = n
    if (nlp>1 .and. nlp< nring) then
      im = n-1; ip = n+1
    elseif (nlp.eq.1) then
      im = n+nring-1; ip = n+1
    elseif (nlp.eq.nring)then
      im = n-1; ip = n-nring+1
    endif
    xq=xb(n); xp1=xb(ip); xm1=xb(im)
    yq=yb(n); yp1=yb(ip); ym1=yb(im)
    !
    !cx=xp1-two*xq+xm1; cx=cx*rds*rds
    !cy=yp1-two*yq+ym1; cy=cy*rds*rds
!
    sp=sqrt((xp1-xq)*(xp1-xq) + (yp1-yq)*(yp1-yq))
    sm=sqrt((xm1-xq)*(xm1-xq) + (ym1-yq)*(ym1-yq))
!
    if (consttension) then
    !constant tension: T=const; tau={dX/dsigma}/|dX/dsigma|; F=d(T*tau)/dsigma
      !version 1: taking derivative for F and then discretize, need cx, cy above, and gd bellow
      !gd=(sp+sm)/(dls*two)
      !Fvx(n)= sw(il)*(cx/gd-half*(xp1-xm1)*rds/(gd*gd)*(sp-sm)*rds)
      !Fvy(n)= sw(il)*(cy/gd-half*(yp1-ym1)*rds/(gd*gd)*(sp-sm)*rds)
      !version 2: discretize F directly
      Fvx(n)= sw(il)*rds*((xp1-xq)/sp-(xq-xm1)/sm)
      Fvy(n)= sw(il)*rds*((yp1-yq)/sp-(yq-ym1)/sm)
    else
    !!
! w rest metric: T=k_e(|dX/dsigma|-J_ref); F=d(T*tau)/dsigma
      Fvx(n)=sw(il)*rds*((one-membrane_reference_metric*dls/sp)*rds* &
           (xp1-xq)-(one-membrane_reference_metric*dls/sm)*rds*(xq-xm1))
      Fvy(n)=sw(il)*rds*((one-membrane_reference_metric*dls/sp)*rds* &
           (yp1-yq)-(one-membrane_reference_metric*dls/sm)*rds*(yq-ym1))
    endif
!=======================================================================
! no bending force
!=======================================================================
  enddo
!
  return
end subroutine IBForceVec
!
end module IBforce
