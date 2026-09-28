!include subroutines for solving Stokes equation, reading parameters and 
! save (u,v,p) data
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscsnes.h>
!
module linsys
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use petscksp
  use petscmat
  use petscvec
  use petscsys

  use parameters
  use IBmod
  use IBforce
  use myfft

  implicit none

  double precision :: myalpha, mybeta, mygam

  public :: myalpha, mybeta, mygam
  public :: velRHS, linmatk, sollinsys, rsoltouvp, asol, inituvp, cmpuvp, &
       readpar, outuvpx, setbc,TransRHS, wrapLinSolve
!
contains
!--------------------------------------------------------------
!
  subroutine readpar(CFL,nu,bke,bki,vco,ntmax,nfreq,dlt &
#ifdef SIMCELL_TESTING
      ,dualchem_diffusion_test_override &
#endif
      )
  integer :: ntmax, nfreq, read_status
  double precision :: CFL,nu,bke,bki,vco,dlt
  real(dp) :: physical_scales(5)
  real(dp) :: expected_velocity_scale,scale_tolerance
  character(len=256) :: read_message
  character(len=256) :: parameter_file
#ifdef SIMCELL_TESTING
  real(dp),intent(in),optional::dualchem_diffusion_test_override
#endif
  namelist/PARAM/CFL,nu,bke,bki,vco,ntmax,nfreq,dlt, &
       dualchem_kp,dualchem_pump_start_time,dualchem_pump_width,dualchem_initial_concentration, &
       dualchem_rear_width,dualchem_rear_amplitude_ratio,stage14_active_actin_feedback, &
       stage12_initial_network_concentration,stage12_initial_free_concentration, &
       enforce_fixed_timestep, &
       dualchem_initial_polarization,dualchem_diffusion, &
       stage12_actin_gamma,stage12_actin_eta,stage12_actin_eta_s, &
       stage12_actin_k_sigma,stage12_actin_dc,stage12_actin_jc, &
       stage12_actin_dw,stage12_actin_theta0,kw,dualchem_kc, &
       stage12_localized_polymerization, &
       stage12_pnas_balanced_actin_profile, &
       stage14_remove_uniform_actin_stress, &
       stage14_recenter_interface,fsi_max_outer_iterations, &
       stage14_shape_force_scale,stage14_adhesion,stage14_external_load_x, &
       water_stress_mobility,water_osmotic_mobility,enable_osmotic_feedback, &
       initial_cell_radius,initial_cell_axis_ratio,initial_cell_shape_factor, &
       initial_cell_center_offset_x,initial_cell_center_offset_y, &
       membrane_reference_metric,clstiff,dualchem_interface_gmres_rtol, &
       scale_length_um,scale_concentration_millimolar,scale_velocity_um_s, &
       scale_time_s,scale_stress_pa
!
  ! Reset before every read so a missing entry cannot inherit a value from a
  ! prior readpar call in a test, restart, or embedding process.
  dualchem_diffusion=zero
  parameter_file='input.par'
  if(command_argument_count()>=1)call get_command_argument(1,parameter_file)
  open(unit=99,file=trim(parameter_file),status='old',form='formatted')
  read(99,nml=PARAM,iostat=read_status,iomsg=read_message)
  close(unit=99)
  if(read_status/=0)then
    write(*,'(a,1x,a)')'readpar: invalid PARAM namelist:',trim(read_message)
    error stop 'readpar: invalid PARAM namelist'
  endif
#ifdef SIMCELL_TESTING
  if(present(dualchem_diffusion_test_override)) &
      dualchem_diffusion=dualchem_diffusion_test_override
#endif
  if(.not.is_finite_run_scalar(dualchem_diffusion))then
    error stop 'readpar: dualchem_diffusion must be finite and positive'
  endif
  if(dualchem_diffusion<=zero)then
    error stop 'readpar: dualchem_diffusion must be finite and positive'
  endif
  if(.not.ieee_is_finite(dualchem_pump_width) .or. dualchem_pump_width<=zero) &
      error stop 'Pump width must be finite and positive'
  if(.not.ieee_is_finite(dualchem_rear_width) .or. dualchem_rear_width==zero) &
      error stop 'Rear width must be positive or negative legacy sentinel'
  if(.not.ieee_is_finite(dualchem_rear_amplitude_ratio) .or. dualchem_rear_amplitude_ratio<zero) &
      error stop 'Rear amplitude ratio must be finite and nonnegative'

  if(.not.ieee_is_finite(stage12_initial_network_concentration) .or. &
       stage12_initial_network_concentration<zero .or. &
       .not.ieee_is_finite(stage12_initial_free_concentration) .or. &
       stage12_initial_free_concentration<zero) &
      error stop 'Initial actin concentrations must be finite and nonnegative'
  used_legacy_water_mobility_fallback=.false.
  if(water_stress_mobility<zero)then
    water_stress_mobility=kw(1)
    used_legacy_water_mobility_fallback=.true.
  endif
  if(water_osmotic_mobility<zero)then
    water_osmotic_mobility=kw(1)
    used_legacy_water_mobility_fallback=.true.
  endif
  if(.not.ieee_is_finite(clstiff) .or. clstiff<zero)then
    error stop 'Dimensionless membrane elastic stiffness must be nonnegative'
  endif
  if(.not.ieee_is_finite(membrane_reference_metric) .or. &
       membrane_reference_metric<zero)then
    error stop 'Membrane reference metric must be finite and nonnegative'
  endif
  if(.not.ieee_is_finite(initial_cell_center_offset_x) .or. &
       .not.ieee_is_finite(initial_cell_center_offset_y) .or. &
       .not.ieee_is_finite(initial_cell_shape_factor))then
    error stop 'Initial cell geometry parameters must be finite'
  endif
  sw(1)=clstiff
  if(.not.ieee_is_finite(stage14_adhesion) .or. stage14_adhesion<zero)then
    error stop 'Dimensionless adhesion A_d must be finite and nonnegative'
  endif
  if(.not.ieee_is_finite(stage14_external_load_x))then
    error stop 'Dimensionless external load must be finite'
  endif
  if(fsi_max_outer_iterations<1)then
    error stop 'FSI outer iteration count must be positive'
  endif
  if(.not.ieee_is_finite(dualchem_interface_gmres_rtol) .or. &
       dualchem_interface_gmres_rtol<=zero .or. &
       dualchem_interface_gmres_rtol>=one)then
    error stop 'Dual-chemical interface GMRES tolerance must lie in (0,1)'
  endif
  if(.not.ieee_is_finite(dualchem_kc) .or. dualchem_kc<zero .or. &
       .not.ieee_is_finite(dualchem_kp))then
    error stop 'Dual-chemical passive and active rates are invalid'
  endif
  if(.not.ieee_is_finite(dualchem_initial_concentration) .or. &
       dualchem_initial_concentration<zero .or. &
       .not.ieee_is_finite(dualchem_initial_polarization))then
    error stop 'Dual-chemical initial data are invalid'
  endif
  if(.not.ieee_is_finite(water_stress_mobility) .or. &
       water_stress_mobility<zero .or. &
       .not.ieee_is_finite(water_osmotic_mobility) .or. &
       water_osmotic_mobility<zero)then
    error stop 'Water mobilities must be finite and nonnegative'
  endif
  physical_scales=[scale_length_um,scale_concentration_millimolar, &
       scale_velocity_um_s,scale_time_s,scale_stress_pa]
  if(.not.all(ieee_is_finite(physical_scales)) .or. &
       any(physical_scales<=zero))then
    error stop 'All physical reference scales must be finite and positive'
  endif
  expected_velocity_scale=scale_length_um/scale_time_s
  scale_tolerance=256.0_dp*epsilon(one)*max(abs(scale_velocity_um_s), &
       abs(expected_velocity_scale),tiny(one))
  if(abs(scale_velocity_um_s-expected_velocity_scale)>scale_tolerance)then
    error stop 'Physical scales must satisfy U0=L0/T0'
  endif
!
  return
  end subroutine readpar

!========================================================================
subroutine wrapLinSolve(f,g,u,v,p, alpha,beta,gam)
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in) :: f, g
  double precision, dimension(-1:nx+1,-1:ny+1), intent(out) :: u,v,p
  double precision :: alpha, beta, gam
!
  complex*16, dimension(-1:nx+1,-1:ny+1) :: rsu,rsv
  double precision, dimension(3*ny-1,nx) :: rsol
  double precision, dimension(-1:nx+1,4) :: homogeneous_bc
!
  call TransRHS(f,g, rsu,rsv) ! transform RHS in (ma1,mb1) to Fourier space
  call sollinsys(rsu,rsv,rsol, alpha, beta, gam) ! real valued solution  in rsol
  call rsoltouvp(rsol,u,v,p)
  ! rsoltouvp fills only the physical staggered unknowns.  Complete every
  ! solve with the boundary representation consumed by interpolation/output:
  ! periodic left/right and strict homogeneous no-slip at south/north.
  homogeneous_bc = zero
  call setbc(u,homogeneous_bc,1)
  call setbc(v,homogeneous_bc,2)
  call setbc(p,homogeneous_bc,3)
!
  return
end subroutine wrapLinSolve
!========================================================================
subroutine setbc(rho,uvbc,isel)
! periodic on left/right, Dirichlet on top/bottom
!
  integer :: isel
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in out) :: rho
!
  integer :: i, j, k
  double precision :: tmp
  double precision, dimension(-1:nx+1) :: top, bot
!
  select case (isel)
  case (1,-1) !u
    if (isel .eq. 1) then
      bot(0:nx-1) = uvbc(0:nx-1,1); top(0:nx-1) = uvbc(0:nx-1,2)
    else ! homogeneous Dirichlet
      bot(0:nx-1) = 0.d0     ; top(0:nx-1) = 0.d0     
    endif
    j = 0
    rho(0:nx-1,j) = two*bot(0:nx-1) - rho(0:nx-1,j+1)
    j = -1
    rho(0:nx-1,j) = two*bot(0:nx-1) - rho(0:nx-1,j+2)
    j = ny+1
    rho(0:nx-1,j) = two*top(0:nx-1) - rho(0:nx-1,j-1)
    !
    i = -1
    rho(i,0:ny+1) = rho(nx+i,0:ny+1)
    i = nx
    rho(i,0:ny+1) = rho(i-nx,0:ny+1)
    i = nx+1
    rho(i,0:ny+1) = rho(i-nx,0:ny+1)
  case (2,-2) !v
    if (isel .eq. 2) then
      bot(1:nx) = uvbc(1:nx,3); top(1:nx) = uvbc(1:nx,4)
    else ! homogeneous Dirichlet
      bot(1:nx) = 0.d0     ; top(1:nx) = 0.d0     
    endif
    j = 0 ! bottom Bdry
    rho(1:nx,j) = bot(1:nx);
    j = -1
    rho(1:nx,j) = two*bot(1:nx) - rho(1:nx,j+2)
    j = ny ! top Bdry
    rho(1:nx,j) = top(1:nx)
    j = ny+1
    rho(1:nx,j) = two*top(1:nx) - rho(1:nx,j-2)
    !
    i = 0
    rho(i,-1:ny+1) = rho(nx+i,-1:ny+1)
!!    i = -1
!!    rho(i,-1:ny+1) = rho(nx+i,-1:ny+1)
    i = nx+1
    rho(i,-1:ny+1) = rho(i-nx,-1:ny+1)
  case (3) !p
    i = 0
    rho(i, 1:ny) = rho(nx+i, 1:ny)
    i = nx+1
    rho(i, 1:ny) = rho(i-nx, 1:ny)
  case default
    print *,'Can only choose (-)1 (u), or (-)2(v)! stop'
    stop
  end select 
!
  return
end subroutine setbc
!========================================================================
subroutine lapop(lap, rho, isel)
! periodic BC on left/right, Homogeneous Dirichlet on top/bottom 
! need advu(0:nx-1,1:ny), advv(1:nx,1:ny-1)
!
  integer :: isel
  double precision, dimension(-1:nx+1,-1:ny+1) :: lap, rho
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
!
  integer :: i, j, istr, iend, jstr, jend
  double precision :: tp1, tp2, tmp, h, dx, dy, rhh
!
  rhh = one/(hg*hg) 
  lap = 0.d0
!
  select case (isel) 
   ! u component                                            
  case (0) 
    istr = 0 
    iend = nx-1 
    jstr = 1 
    jend = ny
   ! v component                                            
  case (1) 
    istr = 1 
    iend = nx
    jstr = 1 
    jend = ny-1 
   ! rho @ cell centers                                     
  case (2) 
    istr = 1 
    iend = nx
    jstr = 1 
    jend = ny
  case default 
    print *, isel 
    print *, 'no such option in lap, stop' 
    stop 
  end select 
  do j = jstr, jend 
  do i = istr, iend 
    lap(i,j) = (rho(i+1,j)+rho(i-1,j)-two*rho(i,j))          &
     &       + (rho(i,j+1)+rho(i,j-1)-two*rho(i,j))          
  enddo 
  enddo 
  lap = rhh*lap
!
  return
end subroutine lapop
!========================================================================
subroutine getAdv(advu, advv, u, v, uvbc)
! periodic BC on left/right, Homogeneous Dirichlet on top/bottom 
! need advu(0:nx-1,1:ny), advv(1:nx,1:ny-1)
!
  double precision, dimension(-1:nx+1,-1:ny+1) :: u, v
  double precision, dimension(-1:nx+1,-1:ny+1), intent(out) :: advu, advv
  double precision, dimension(-1:nx+1,4) :: uvbc !ubt, utp, vbt, vtp
!
  integer :: i, j, k
  double precision :: tp1, tp2, tmp, h, dx, dy
  double precision, dimension(-1:nx+1,-1:ny+1) :: ua, va
!
  h = hg; dx = h; dy = h
  advu = 0.d0; advv = 0.d0;
! 
!!  call setbc(u,uvbc,1) ! ghost of u CHANGED here
!!  call setbc(v,uvbc,2) ! ghost of v CHANGED here
!
  do i = 0,nx-1 ! average of v @ u locations
    do j = 1, ny
      va(i,j) = 0.25d0*(v(i,j-1) + v(i+1,j-1) + v(i,j) + v(i+1,j))
    enddo
  enddo
  call setbc(va,uvbc,1) ! need to be modified to account for the shifted BC loc
  do i = 1, nx ! average of u @ v locations
    do j = 1, ny-1
      ua(i,j) = 0.25d0*(u(i-1,j) + u(i,j) + u(i-1,j+1) + u(i,j+1))
    enddo
  enddo
  call setbc(ua,uvbc,2) ! 
!
  do i = 0, nx-1
    do j = 1, ny
      advu(i,j) = u(i,j)*(u(i+1,j)-u(i-1,j))+va(i,j)*(u(i,j+1)-u(i,j-1)) + & 
                  (u(i+1,j)*u(i+1,j)  - u(i-1,j)*u(i-1,j) ) +  &
                  (u(i,j+1)*va(i,j+1) - u(i,j-1)*va(i,j-1))
    enddo
  enddo
  advu = advu*half*half/h
  do i = 1, nx
    do j = 1, ny-1
      advv(i,j) = ua(i,j)*(v(i+1,j)-v(i-1,j))+v(i,j)*(v(i,j+1)-v(i,j-1)) + &
                  (ua(i+1,j)*v(i+1,j) - ua(i-1,j)*v(i-1,j)) +  &
                  (v(i,j+1)*v(i,j+1)  - v(i,j-1)*v(i,j-1))
    enddo
  enddo
  advv = advv*half*half/h
!
  return
end subroutine getAdv
!========================================================================
!========================================================================
subroutine TransRHS(f,g, rsu,rsv)
! Transform RHS vars (f,g) to complex space (rsu,rsv)
  complex*16, dimension(-1:nx+1,-1:ny+1), intent(out) :: rsu, rsv
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in) :: f, g
  !
  integer :: j, nl
  complex*16, dimension(1:nx) :: din, dout
  integer*8 :: plan
!
  nl = nx
!
  do j = 1, ny
    din = f(0:nx-1,j)
    call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,din,dout)
    rsu(0:nx-1,j) = dout(1:nx)
    call dfftw_destroy_plan(plan)
  enddo
!========================================================================
!========================================================================
  ! FFT along each row (j=1:ny-1) for v-component
  do j = 1, ny-1
    din = g(1:nx,j)
    call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,din,dout)
    rsv(1:nx,j) = dout(1:nx)
    call dfftw_destroy_plan(plan)
  enddo
!
  return
end subroutine TransRHS
!========================================================================
!========================================================================
subroutine velRHS(u, v, xb,yb, f, g)
! get RHS for Stokes solver, no transformation performed
! no sigma_n here, will need to add for two phase case
  implicit none
!
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in) :: u, v
  double precision, dimension(mpts), intent(in) :: xb,yb
  double precision, dimension(-1:nx+1,-1:ny+1), intent(out) :: f, g
!
  double precision, dimension(7,mpts) :: mpoly
  double precision, dimension(mpts) :: fsx, fsy, nvx,nvy
  double precision, dimension(-1:nx+1,-1:ny+1):: tiu, tiv, teu, tev
!
  !jdu & jdv need to be checked for current (xb,yb)!!
!  call getNormal(xb,yb, nvx,nvy, mpoly) !normal (nvx,nvy) needed for sgmn
  !call IBJacobian(xb,yb, jv)

  tiu = bki*dble(jdu(-1:nx+1,-1:ny+1))*u
  tiv = bki*dble(jdv(-1:nx+1,-1:ny+1))*v
  teu = bke*dble(1-jdu(-1:nx+1,-1:ny+1))*u
  tev = bke*dble(1-jdv(-1:nx+1,-1:ny+1))*v
!
  call IBForceVec(xb,yb, fsx,fsy)
  call newSpread(xb,yb,fsx,fsy,f,g)
  f = f-(tiu+teu)
  g = g-(tiv+tev)

!noTrans  do j = 1, ny
!noTrans    din = advu(0:nx-1,j)
!noTrans    call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
!noTrans    call dfftw_execute_dft(plan,din,dout)
!noTrans    rsu(0:nx-1,j) = dout(1:nx)
!noTrans    call dfftw_destroy_plan(plan)
!noTrans  enddo
!========================================================================
! modify the RHS due to Dirichlet BC from u-component @ bottom
!!  din = uvbc(0:nx-1,1)
!!  call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
!!  call dfftw_execute_dft(plan,din,dout)
!!  rsu(0:nx-1,1) = rsu(0:nx-1,1)+two*beta*dout(1:nx)
!!  call dfftw_destroy_plan(plan)
!!! modify the RHS due to Dirichlet BC from u-component @ top
!!  din = uvbc(0:nx-1,2)
!!  call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
!!  call dfftw_execute_dft(plan,din,dout)
!!  rsu(0:nx-1,ny) = rsu(0:nx-1,ny)+two*beta*dout(1:nx)
!!  call dfftw_destroy_plan(plan)
!========================================================================
  ! FFT along each row (j=1:ny-1) for v-component
!noTrans  do j = 1, ny-1
!noTrans    !!din = 0.d0; 
!noTrans    din = advv(1:nx,j)
!noTrans    call dfftw_plan_dft_1d(plan, nl, din,dout,FFTW_FORWARD,FFTW_ESTIMATE)
!noTrans    call dfftw_execute_dft(plan,din,dout)
!noTrans    rsv(1:nx,j) = dout(1:nx)
!noTrans    call dfftw_destroy_plan(plan)
!noTrans  enddo
!
  return
end subroutine velRHS

!--------------------------------------------------------------         
!--------------------------------------------------------------         
subroutine linmatk(AB, k, alpha, beta, gam)
  implicit none
  !!include '/Users/yaol/local/fftw/include/fftw3.f'
  integer :: isel, k ! k is the wave number
  double precision :: alpha, beta, gam
  ! LDAB = 2*3+3+1 = 10
  complex*16, dimension(10,3*ny-1) :: AB
!
  integer :: i, j, jj, ik, jk, kl, ku, nl
  double precision :: tpikh, pikh, cstpikh, sitpikh, cspikh, sipikh
  double precision :: h, alpbet, alpbet1
  complex*16 :: Epikh, Etpikh, Etpikh1, Etpikhp
  complex*16, dimension(3) :: l1, l2, l3, u1, u2, u3
  complex*16, dimension(7) :: c1, c2, c3
!
  nl = 3*ny-1;  h = hg
  ! The FFT mode is indexed on nx periodic grid points, so its one-cell phase
  ! is theta_k=2*pi*k/nx.  Multiplying by the physical spacing hg is only
  ! equivalent on a unit-length domain; it doubled every phase after moving
  ! the shared domain to [-1,1]^2 and corrupted the Stokes symbol.
  tpikh = two*cpi*dble(k)/dble(nx)
  cstpikh = dcos(tpikh); sitpikh = dsin(tpikh)
!
  alpbet = alpha-2.d0*beta*cstpikh
  alpbet1= alpha-2.d0*beta*cstpikh+beta
!
  ! CMPLX without KIND returns default complex even when its arguments are
  ! double precision.  Preserve the COMPLEX*16 Stokes symbol coefficients;
  ! the single-precision round trip is large enough to raise the FSI residual
  ! floor on the two-unit domain.
  Etpikh1= cmplx(cstpikh-1.d0,sitpikh,kind=dp)*gam
  Etpikhp= cmplx(1.d0-cstpikh,sitpikh,kind=dp)*gam
!
  AB = cmplx(zero,zero,kind=dp)
!
  u1(1) = 0.d0;     u1(2) = 0.d0;  u1(3) =-gam
  l1(1) = Etpikhp;  l1(2) = gam;   l1(3) = 0.d0
  u2(1) =-beta;     u2(2) = 0.d0;  u2(3) = Etpikh1
  l2(1) = 0.d0;     l2(2) = 0.d0;  l2(3) =-beta
  u3(1) =-beta;     u3(2) =-gam  ; u3(3) = 0.d0
  l3(1) = gam  ;    l3(2) = 0.d0;  l3(3) =-beta
  c1(1:3) = u1; c1(4) = cmplx(0.d0,0.d0,kind=dp);   c1(5:7) = l1
  c2(1:3) = u2; c2(4) = cmplx(alpbet,0.d0,kind=dp); c2(5:7) = l2
  c3(1:3) = u3; c3(4) = cmplx(alpbet,0.d0,kind=dp); c3(5:7) = l3
!
  kl = 3; ku = 3
!
  jj = 1;
  j  = 3*(jj-1)+1; i = 1; ik = kl + ku + 1 + i -j
  AB(ik:ik+3,j) = c1(4:7)
  if (k .eq. 0) then
    AB(ik,j) = cmplx(1.0d0,0.d0,kind=dp)
  endif
  j  = 3*(jj-1)+2; i = 1; ik = kl + ku + 1 + i -j
  AB(ik:ik+4,j) = c2(3:7)
  AB(ik+1,j) = cmplx(alpbet1,0.d0,kind=dp)
  j  = 3*(jj-1)+3; i = 1; ik = kl + ku + 1 + i -j
  AB(ik:ik+5,j)= c3(2:7)
  do jj = 2, ny-2
    j = 3*(jj-1)+1; i = j-3; ik = kl + ku + 1 + i -j
    AB(ik:ik+6,j) = c1
    j = 3*(jj-1)+2; i = j-3; ik = kl + ku + 1 + i -j
    AB(ik:ik+6,j) = c2
    j = 3*(jj-1)+3; i = j-3; ik = kl + ku + 1 + i -j
    AB(ik:ik+6,j) = c3
  enddo
  jj= ny-1
  j = 3*(jj-1)+1; i = j-3; ik = kl + ku + 1 + i -j
  AB(ik:ik+6,j) = c1
  j = 3*(jj-1)+2; i = j-3; ik = kl + ku + 1 + i -j
  AB(ik:ik+6,j) = c2
  j = 3*(jj-1)+3; i = j-3; ik = kl + ku + 1 + i -j
  AB(ik:ik+5,j) = c3(1:6)
!
  jj= ny
  j = 3*(jj-1)+1; i = j-3; ik = kl + ku + 1 + i -j 
  AB(ik:ik+4,j) = c1(1:5)
  j = 3*(jj-1)+2; i = j-3; ik = kl + ku + 1 + i -j 
  AB(ik:ik+2,j) = c2(1:3)
  AB(ik+3,j) = cmplx(alpbet1,0.d0,kind=dp)
!
  !!prt if (k .eq. 0) then
  !!prt print *, 'Coefficient Matrix start'
  !!prt do i = 1, nl
  !!prt   !!print '(i3, 1024("(",e8.1,1x,e8.1,")"))', i, AB(4:10,i)
  !!prt   print '(i3, 1024(e9.2,1x))', i, real(AB(4:10,i))
  !!prt enddo
  !!prt print *, 'Coefficient Matrix end'
  !!prt endif
!
  return
end subroutine linmatk
!--------------------------------------------------------------         
subroutine rsoltouvp(rsol, u, v, p)
  double precision, dimension(3*ny-1,nx) :: rsol
  double precision, dimension(-1:nx+1,-1:ny+1) :: u,v,p
  integer ::j, jj, k
!
!!  u = 0.; v= 0.; p = 0.
!
  !!do k = 1, nx
  !!  do jj = 1, ny-1
  !!    j = 3*(jj-1)+1
  !!    p(k,  ny-jj+1) = rsol(j,  k)
  !!    u(k-1,ny-jj+1) = rsol(j+1,k)
  !!    v(k,  ny-jj) = rsol(j+2,k)
  !!  enddo
  !!  jj = ny
  !!  j = 3*(jj-1)+1
  !!  p(k,  ny-jj+1) = rsol(j,  k)
  !!  u(k-1,ny-jj+1) = rsol(j+1,k)
  !!enddo
  do k = 1, nx
    do jj = 1, ny-1
      j = 3*(jj-1)+1
      p(k,  ny-jj+1) = rsol(j,  k)
      u(k-1,ny-jj+1) = rsol(j+1,k)
      v(k,  ny-jj) = rsol(j+2,k)
    enddo
    jj = ny
    j = 3*(jj-1)+1
    p(k,  ny-jj+1) = rsol(j,  k)
    u(k-1,ny-jj+1) = rsol(j+1,k)
  enddo
!
  return
end subroutine rsoltouvp


subroutine sollinsys(rsu,rsv,rsol, alpha, beta, gam)
  implicit none
  complex*16, dimension(-1:nx+1,-1:ny+1), intent(in) :: rsu, rsv
  double precision, dimension(3*ny-1,nx), intent(out) :: rsol
!  double precision :: dt
!
  integer :: i, j, k, iflg, info, nl
!  double precision :: dx, tupih, pih!, h
  double precision :: alpha, beta, gam
  complex*16, dimension(10,3*ny-1) :: AB, AFB
  complex*16, dimension(3*ny-1) :: rhs, solv, solw
  complex*16, dimension(3*ny-1,nx) :: resul!, csol!!, rhsm
  integer :: ipiv(3*ny-1), nlen
  character :: equed
  double precision :: mvR(3*ny-1), mvC(3*ny-1), ferr(1), berr(1), rwork(6*ny-2)
  double precision :: rcond
  complex*16 :: work(6*ny-2)
  integer*8 :: plan
  complex*16, dimension(nx) :: din, dout

!
!!  h = hg
!!prt print *, '====='
!!prt   print*,alpha, beta, gam
!!prt print *, '====='
!
  nl = 3*ny-1
!
  do k = nx, 1, -1 !need to have 2:nx before compute k=1
  !!do k = 1, nx
    rhs = cmplx(zero,zero,kind=dp)
    do i = 1, ny-1
      j = 3*(i-1) + 1
      rhs(j  ) = cmplx(0.d0,0.d0,kind=dp)
      rhs(j+1) = rsu(k-1,ny-i+1)
      rhs(j+2) = rsv(k,  ny-i)
    enddo
    i = ny
    j = 3*(i-1)+1
    rhs(j) = cmplx(0.d0,0.d0,kind=dp)
    rhs(j+1) = rsu(k-1,ny-i+1)
!!      print '(1x,"about??? ", i5, 1024(e14.6,1x))', k, maxval(real(rhs)), &
!!      minval(real(rhs))
!
!!    rhsm(:,k) = rhs
!
    AB = zero
!
    !!call linmatk(AB, k-nx/2+1, alpha, beta, gam, 0)
    call linmatk(AB, k-1, alpha, beta, gam)
!!    print '(1x,"howabout AB ", 1024(e14.6,1x))', maxval(real(AB)), &
!!    minval(real(AB))
    info = 0
    ipiv = 0

    solw = rhs
!    call zgbsv(nl, 3, 3, 1, AB, 10, ipiv, rhs, nl, info)
            if (k .eq. 1) then ! for singular matrix at 0 wave number
              call zgbsv(nl-1, 3, 3, 1, AB(:,2:nl), 10, ipiv, rhs(2:nl), nl-1, info)
              rhs(1) = -sum(resul(1,2:nx))
            else
              call zgbsv(nl, 3, 3, 1, AB, 10, ipiv, rhs, nl, info)
            endif
!!deb    AFB = AB
!!deb    equed='N'
!!deb    call ZGBTRF(nl, nl, 3, 3, AFB, 10, ipiv, info )
!!deb!!DEBUG
!!deb    IF( INFO.NE.0 ) THEN
!!deb      print '(1x,1024(e14.6,1x))', AFB
!!deb      stop
!!deb    ENDIF
!!deb!!DEBUG
!!deb    mvR = zero; mvC = zero; rhs = solw
!!deb    call zgbsvx('F','N',nl,3,3,1,AB,10,AFB,nl,ipiv,&
!!deb       equed,mvR,mvC,rhs,nl,solv,nl,rcond,ferr,berr,work,rwork,info)
!!deb    rhs = solv
!!    if (k .eq. 2) then
!!      do i = 1, 3*ny-1
!!      print *, rhs(i)
!!      enddo
!!    endif
    if (info .eq. 0) then
      resul(1:nl,k) = rhs
!!      print '(1x,"get rhs? ", i5, 1024(e14.6,1x))', k, maxval(real(rhs)), &
!!      minval(real(rhs)), maxval(dimag(rhs)), maxval(dimag(rhs))
    else
      print *, nl
      print '("Error in wave stop ", 2(i5,1x) 10(e22.15,1x))', k, info, rcond
!!      print'(1x, 1024("(",e11.3,1x,e11.3,")"))', sum(AB(4:10,:),1)
      print'(1x, 1024(e11.3,1x))', (sum(AB(4:10,:),1))
      stop
    endif
  enddo
  rhs = zero
!

  rsol = 0.d0; !csol = 0.d0
  nlen = nx
  do i = 1, 3*ny-1
    dout = resul(i,1:nx); 
!!    print '("dout ", i5,1024(e14.6,1x))',i, maxval(real(dout)), minval(real(dout))
    call dfftw_plan_dft_1d(plan, nlen, dout,din,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,dout,din)
    !csol(i,:) = din/dble(nlen)
    rsol(i,:) = real(din/dble(nlen))
!!    print '("din  ", i5,1024(e14.6,1x))',i, maxval(real(din)), minval(real(din)),&
!!       maxval(dimag(csol)), minval(dimag(csol))
    call dfftw_destroy_plan(plan)
  enddo
!!prt  print *, norm2(real(csol)), norm2(imag(csol))
!!prt  do i = 1, 3*ny-1
!!prt  !!do i = 1, 2
!!prt    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, resul(i,:)
!!prt    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, csol(i,:)
!!prt    print'(1x, 1024(e11.3,1x))', rsol(i,:)
!!prt  enddo
!!  print *, 'rhs'
!!  do i = 1, 3*ny-1
!!  !!do i = 1, 2
!!    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, resul(i,:)
!!    print'(1x, i3, 1024(e11.3,1x))', i, real(resul(i,:))
!!    !!print'(1x, i3, 1024("(",e11.3,1x,e11.3,")"))', i, rhsm(i,:)
!!    !!print'(1x, i3, 1024(e11.3,1x))', i, real(rhsm(i,:))
!!  enddo
!
  return
end subroutine sollinsys
!--------------------------------------------------------------         
!!double precision function asol(x,y,t,isel)
subroutine asol(mysol, x,y,t,isel)
!!  use parameters
  implicit none
!
  double precision :: x, y, t, mysol
  integer :: isel
!
  tupi = two*cpi; pi = cpi
  select case (isel) ! first two rows are for steady Stokes
  case (0) !u
    mysol =-sin(tupi*y)*sin(pi*x)*sin(pi*x)*exp(-t)
    !!case4 mysol = cos(tupi*x)*(three*y*y-two*y)*exp(-t)
    !!case2 mysol =-sin(tupi*y)*sin(pi*x)*sin(pi*x)*exp(-t)
    !!case1 mysol = cos(tupi*x)*(three*y*y-two*y)*exp(-t)
    !!case3 mysol = 0.d0
  case (1) !v
    mysol = sin(tupi*x)*sin(pi*y)*sin(pi*y)*exp(-t)
    !!case4 mysol = tupi*sin(tupi*x)*y*y*(y-one)*exp(-t)
    !!case2 mysol = sin(tupi*x)*sin(pi*y)*sin(pi*y)*exp(-t)
    !!case1 mysol = tupi*sin(tupi*x)*y*y*(y-one)*exp(-t)
    !!case3 mysol = 0.d0
!!  case (-2) !px
!!    !!mysol =cos(tupi*x)*cos(tupi*y)
!!    mysol = tupi*pi*cos(tupi*x)*sin(tupi*y)
!!    !!mysol = 0.d0
!!  case (-3) !py
!!    !!mysol =cos(tupi*x)*cos(tupi*y)
!!    mysol = tupi*pi*sin(tupi*x)*cos(tupi*y)
!!    !!mysol = 0.d0
  case (2) !p
    mysol = pi*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case4 mysol = nu*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case2 mysol = pi*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case1 mysol = nu*sin(tupi*x)*sin(tupi*y)*exp(-t)
    !!case3 mysol = 0.d0
  case (3) ! forcing on x
    !!mysol = two*pi*pi*(-nu+( one +two*nu)*cos(tupi*x))*sin(tupi*y)*exp(-t)
    mysol = half*(one-four*pi*pi*nu+(-one +pi*pi*four*( one+two*nu))*cos(tupi*x))*sin(tupi*y)*exp(-t)
    !!case4 mysol = two*nu*cos(tupi*x)*(-three+two*pi*pi*y*(-two +three*y)+pi*sin(tupi*y))*exp(-t)
    !!case2 mysol = half*(one-four*pi*pi*nu+(-one +pi*pi*four*( one+two*nu))*cos(tupi*x))*sin(tupi*y)*exp(-t)
    !!case1 mysol = exp(-t)*cos(tupi*x)*(-6.d0*nu+y*two*(one-four*pi*pi*nu)+three*y*y*(-one+tupi*tupi*nu)+tupi*nu*sin(tupi*y))
    !!case3 mysol = 0.d0
  case (4) ! forcing on y
    !!mysol =-two*pi*pi*(-nu+(-one +two*nu)*cos(tupi*y))*sin(tupi*x)*exp(-t)
    mysol =-half*(one-four*pi*pi*nu+(-one +pi*pi*four*(-one+two*nu))*cos(tupi*y))*sin(tupi*x)*exp(-t)
    !!case4 mysol = two*pi*nu*(two-6.d0*y-four*pi*pi*y*y+four*pi*pi*y*y*y+cos(tupi*y))*sin(tupi*x)*exp(-t)
    !!case2 mysol =-half*(one-four*pi*pi*nu+(-one +pi*pi*four*(-one+two*nu))*cos(tupi*y))*sin(tupi*x)*exp(-t)
    !!case1 mysol = exp(-t)*tupi*sin(tupi*x)*(two*nu-6.d0*y*nu+y*y*(one-four*pi*pi*nu)+y*y*y*(-one+tupi*tupi*nu)+nu*cos(tupi*y))
    !!case3 mysol = 0.d0
  case default
    print *, 'no this thing'
    stop
  end select
!
  return
END subroutine asol
!
subroutine inituvp(u0, v0, p0, f, g, uvbc, isel, dt, time)
  use parameters
!!  use linsys
  !use linsys
  implicit none
!
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in out) :: u0, v0, p0, f, g
  double precision, dimension(-1:nx+1,4), intent(in out) :: uvbc !ubt, utp, vbt, vtp
  integer :: isel
  double precision :: time, dt
!
  integer :: i, j
  double precision :: x, y, dx, dy, realt
!!  double precision asol
!!  external asol
!
  dx = hg; dy = hg
  tupi = two*cpi; pi = cpi
!
!!  if (isel .eq. 0) then
  do j = 1, ny
    y = (dble(j)-0.5d0)*dy
    do i = 0, nx-1
      x = (dble(i))*dx
      !!u0(i,j) = asol(x,y,time,0)
      call asol(u0(i,j),x,y,time,0)
    enddo
  enddo
  do i = 0,nx-1
    x = (dble(i))*dx
    y = ymin
    !!uvbc(i,1) = asol(x,y,time,0)
    call asol(uvbc(i,1),x,y,time,0)
    y = ymax
    !!uvbc(i,2) = asol(x,y,time,0)
    call asol(uvbc(i,2),x,y,time,0)
  enddo
  do j = 1, ny-1
    y = (dble(j))*dy
    do i = 1, nx
      x = (dble(i)-0.5)*dx
      !!v0(i,j) = asol(x,y,time,1)
      call asol(v0(i,j),x,y,time,1)
    enddo
  enddo
  do i = 1, nx
    x = (dble(i)-0.5)*dx
    y = ymin
    call asol(uvbc(i,3),x,y,time,1)
    y = ymax
    call asol(uvbc(i,4),x,y,time,1)
  enddo
  do j = 1, ny
    y = (dble(j)-.5)*dy
    do i = 1, nx
      x = (dble(i)-0.5)*dx
      !!p0(i,j) = asol(x,y,time,2)
      call asol(p0(i,j),x,y,time,2)
    enddo
  enddo
!!  endif
  select case (isel)
  case (0)
    realt = time
  case (1)
    realt = time - 0.5*dt
!!    realt = time
  end select

    do j = 1, ny
      y = (dble(j)-0.5)*dy
      do i = 0, nx-1
        x = (dble(i)) *dx
        !!f(i,j) = asol(x,y,realt,3)
        call asol(f(i,j),x,y,realt,3)
      enddo
    enddo
    do j = 1, ny-1
      y = (dble(j))*dy
      do i = 1, nx
        x = (dble(i)-0.5)*dx
        !!g(i,j) = asol(x,y,realt,4)
        call asol(g(i,j),x,y,realt,4)
      enddo
    enddo
!!  endif
!!  call setbc(v0, p0, p0, time, isel) 
!
  return
end subroutine inituvp
!
subroutine cmpuvp(u,v,p,time)
  implicit none
!
  double precision, dimension(-1:nx+1,-1:ny+1) :: u,v,p,px, py
  double precision, dimension(-1:nx+1,-1:ny+1) :: ru,rv,rp, rpx, rpy, rqx, rqy
  double precision :: time
  !!double precision sol
  !!external sol
!
  double precision :: eu, ev, ep, dx, dy, x, y, mu, mv, mp, mpx, mpy, epx, epy
  integer :: i, j, k, isel
!
  dx = hg; dy = hg
!
  do j = 1, ny
    y = (dble(j)-0.5)*dy
    do i = 0, nx-1
      x = (dble(i)) *dx
      !!ru(i,j) = sol(x,y,time,0)
      call asol(ru(i,j),x,y,time,0)
    enddo
  enddo
  do j = 1, ny
    y = (dble(j))*dy
    do i = 1, nx
      x = (dble(i)-0.5)*dx
      !!rv(i,j) = sol(x,y,time,1)
      call asol(rv(i,j),x,y,time,1)
    enddo
  enddo
  do j = 1, ny
    y = dble(j-0.5)*dy
    do i = 1, nx
      x = dble(i-0.5)*dx
      !!rp(i,j) = sol(x,y,time,2)
      call asol(rp(i,j),x,y,time,2)
    enddo
  enddo
  eu = norm2(     ru(0:nx-1,1:ny)-u(0:nx-1,1:ny))
  mu = maxval(abs(ru(0:nx-1,1:ny)-u(0:nx-1,1:ny)))
  ev = norm2(     rv(1:nx,1:ny-1)-v(1:nx,1:ny-1))
  mv = maxval(abs(rv(1:nx,1:ny-1)-v(1:nx,1:ny-1)))
  ep = norm2(rp(1:nx,1:ny)-p(1:nx,1:ny))
  mp = maxval(abs(rp(1:nx,1:ny)-p(1:nx,1:ny)))

!!deb  rqx = 0.
!!deb  do j = 1, ny
!!deb    do i = 1, nx-1
!!deb      rqx(i,j) = p(i+1,j)-p(i,j)
!!deb    enddo
!!deb    i = 0
!!deb    rqx(i,j) = p(i+1,j)-p(nx,j)
!!deb  enddo
!!deb  rqx = rqx/dx
!!deb  do j = 1, ny
!!deb    y = (dble(j)-.5)*dy
!!deb    do i = 0, nx-1
!!deb      x = dble(i)*dx
!!deb      rpx(i,j) = sol(x,y,time,-2)
!!deb    enddo
!!deb  enddo
!!deb  rqy = 0.
!!deb  do j = 1, ny-1
!!deb    do i = 1, nx
!!deb      rqy(i,j) = p(i,j+1) - p(i,j)
!!deb    enddo
!!deb  enddo
!!deb  rqy = rqy/dy
!!deb  do j = 1, ny
!!deb    y = (dble(j))*dy
!!deb    do i = 1, nx
!!deb      x = (dble(i)-0.5)*dx
!!deb      rpy(i,j) = sol(x,y,time,-3)
!!deb    enddo
!!deb  enddo
!!deb  epx= norm2(rpx(0:nx-1,1:ny)-rqx(0:nx-1,1:ny))
!!deb  mpx= maxval(abs(rpx(0:nx-1,1:ny)-rqx(0:nx-1,1:ny)))
!!deb  epy= norm2(rpy(1:nx,1:ny-1)-rqy(1:nx,1:ny-1))
!!deb  mpy= maxval(abs(rpy(1:nx,1:ny-1)-rqy(1:nx,1:ny-1)))

  print '("Time@", e11.3, " e(L_2)u,v,p ",3(e14.6,1x), " e(L_inf) ", 3(e14.6,1x))', &
    & time, eu*hg, ev*hg, ep*hg, mu, mv, mp
  !!print '(1x,"Time @", e11.3, " error",1024(e14.6,1x))', time, eu*hg, ev*hg, ep*hg, &
  !!  & epx*hg, epx*hg, mu, mv, mp, mpx, mpy
  !!prt print *, '===================='
!!  print *, ''
!!!!  print *, 'Uc=['
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', u(0:nx-1,i)
!!  enddo
!!!!  print *, ']'
!!!!  print *, 'Vc=['
!!  do i = 1,ny-1
!!    print '(1x,1024(e14.6,1x))', v(1:nx,i)
!!  enddo
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', p(1:nx,i)
!!  enddo
!!!!  print *, ']'
!!!!  print *, 'U=['
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', ru(0:nx-1,i)
!!  enddo
!!!!  print *, ']'
!!!!  print *, 'V=['
!!  do i = 1,ny-1
!!    print '(1x,1024(e14.6,1x))', rv(1:nx,i)
!!  enddo
!!  do i = 1,ny
!!    print '(1x,1024(e14.6,1x))', rp(1:nx,i)
!!  enddo
!!  print *, ']'

!
  return
end subroutine cmpuvp
!========================================================================
subroutine outuvpx(n,u,v,p,xpt,ypt,nct)
  integer :: n, nct
  double precision, dimension(-1:nx+1,-1:ny+1), intent(in) :: u,v,p
  double precision, dimension(n), intent(in) :: xpt, ypt
!
  integer :: i, strlen
  character(40) efile
!
  strlen = len_trim(runname)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.ib.',nct
  open(67,file=efile,form='formatted',action='write')
  do i = 1, n
    write(67,'(2(e22.14,1x))')xpt(i),ypt(i)
  enddo
  close(67)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.u.',nct
  open(67,file=efile,access='stream',action='write')
  write(67)u(0:nx-1,1:ny)
  close(67)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.v.',nct
  open(67,file=efile,access='stream',action='write')
  write(67)v(1:nx,1:ny-1)
  close(67)
  write(efile,'(2a,i4.4)') runname(1:strlen),'.p.',nct
  open(67,file=efile,access='stream',action='write')
  write(67)p(1:nx,1:ny)
  close(67)
! 
  return
end subroutine outuvpx
!

end module linsys
