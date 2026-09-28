module geometry_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use parameters, only: dp, pi, tupi, hg, zero, one, two, half, nx, ny, npts, &
      nib, nfil, xshift, yshift, xmin, ymin, mu, cpi, scc
  use chemical_interface_flux_mod,only: CHEMICAL_FLUX_INVALID,CHEMICAL_FLUX_OK
  use grid_types
  use fsi_dualchem_velocity_bridge_mod, only: VELOCITY_BRIDGE_OK, &
      interpolate_fsi_velocity_to_marker
  use interface_side_mod, only: INTERFACE_SIDE_OK, &
      classify_nearest_normal_point
  use small_solver_mod
  implicit none
  private
  public :: update_geometry_wrapper, getCorrection, compute_correction_on_grid, valatIBpt, evalCx, evalmpoly, apply_fresh_cleared_correction, apply_full_correction
  ! Stage 10 exposes the already existing six-row local solve as a narrow test
  ! and reuse seam.  Its algebra remains unit-diffusion; actin passes the
  ! physically equivalent normalized data v/D, r/D, dt_effective=D*dt, and
  ! source_jump/D.  No legacy chemical call is changed by this export.
  public :: solve2DCauchy
  public :: interp_field_to_lagpoint, interp_field_bilinear
  public :: report_fixed_stencil_audit
  integer(int64),save::fixed_stencil_calls=0_int64
  integer(int64),save::fixed_stencil_rank_failures=0_int64
  integer(int64),save::fixed_stencil_near_rank_deficient=0_int64
  integer(int64),save::fixed_stencil_reconstruction_failures=0_int64
  real(dp),save::fixed_stencil_min_rcond=huge(one)
  real(dp),save::fixed_stencil_max_l1_residual=zero
  real(dp),save::fixed_stencil_max_scaled_l1_residual=zero
#ifdef SIMCELL_TESTING
  public :: get_last_correction_marker_velocity
  real(dp),save::last_correction_marker_velocity(npts,2)=zero
  logical,save::last_correction_marker_velocity_ready=.false.
#endif

contains

  subroutine report_fixed_stencil_audit()
    real(dp)::reported_minimum

    reported_minimum=fixed_stencil_min_rcond
    if(fixed_stencil_calls==0_int64 .or. &
         reported_minimum==huge(one))reported_minimum=zero
    write(*,'(a,1x,4(i0,1x),3(es24.16,1x))') &
         'FIG2_FIXED_STENCIL_AUDIT',fixed_stencil_calls, &
         fixed_stencil_rank_failures,fixed_stencil_near_rank_deficient, &
         fixed_stencil_reconstruction_failures,reported_minimum, &
         fixed_stencil_max_l1_residual, &
         fixed_stencil_max_scaled_l1_residual
  end subroutine report_fixed_stencil_audit

  ! Wrapper to be called by LagrangianGrid%update_geometry or directly
  subroutine update_geometry_wrapper(lag_grid, eul_grid, dt, isel)
    type(LagrangianGrid), intent(inout) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    real(dp), intent(in) :: dt
    integer, intent(in) :: isel
    
    ! 1. Compute smooth representation (BdyQuadParametric)
    call compute_boundary_poly(lag_grid)
    
    ! 2. Tagging
    call tag_grid(lag_grid, eul_grid, isel)
    
    ! 3. Compute geometric links (part of getGeometry)
    call compute_grid_links(lag_grid, eul_grid, isel)
    
  end subroutine update_geometry_wrapper

  subroutine compute_boundary_poly(lag_grid)
    type(LagrangianGrid), intent(inout) :: lag_grid
    
    integer :: k, kprev, knext, idx, idx_next, idx_prev, il
    real(dp) :: s(lag_grid%npts + 1), s0, s1, s2, ds
    real(dp) :: x0, x1, x2, y0, y1, y2
    real(dp) :: tp1, tp2, tmp
    real(dp) :: mA(3,3), AF(3,3), rhs(3), sol(3), mR(3), mC(3), work(12)
    real(dp) :: rcond, ferr(1), berr(1), pc(3), rx(3), ry(3), pA(3,3), xx(3)
    integer :: ipiv(3), info, iwork(3)
    character(1) :: equed
    
    ! External LAPACK routine
    ! External LAPACK routine removed
    
    s(1:lag_grid%npts) = scc(1:lag_grid%npts)
    s(lag_grid%npts+1) = tupi + s(1)
    
    ! Loop over points
    do k = 1, lag_grid%npts
        kprev = merge(lag_grid%npts, k-1, k == 1)
        knext = merge(1, k+1, k == lag_grid%npts)
        
        s0 = s(kprev)
        s1 = s(k)
        s2 = s(knext)
        
        if (k == 1) then
            s0 = s0 - tupi
        else if (k == lag_grid%npts) then
            s2 = tupi
        endif
        
        x0 = lag_grid%x(knext); x1 = lag_grid%x(k); x2 = lag_grid%x(kprev)
        y0 = lag_grid%y(knext); y1 = lag_grid%y(k); y2 = lag_grid%y(kprev)
        
        lag_grid%mk(1, k) = s(k)
        
        ! Setup matrix for x coefficients
        mA(1,1)=1.0_dp; mA(1,2)=s0; mA(1,3)=s0*s0
        mA(2,1)=1.0_dp; mA(2,2)=s1; mA(2,3)=s1*s1
        mA(3,1)=1.0_dp; mA(3,2)=s2; mA(3,3)=s2*s2
        
        rx(1)=x2; rx(2)=x1; rx(3)=x0
        
        ! Solve for x coefficients using solve_linear_system
        pc = rx
        pA = mA
        call solve_linear_system(3, pA, pc, xx, rcond, info)
        
        if (info /= 0) then
            print *, "solve_linear_system failed in compute_boundary_poly (x) at k=", k
            stop
        endif
        
        lag_grid%mk(2:4, k) = xx(1:3)
        
        ! Solve for y coefficients
        ry(1)=y2; ry(2)=y1; ry(3)=y0
        pc = ry
        pA = mA
        call solve_linear_system(3, pA, pc, xx, rcond, info)
                    
        if (info /= 0) then
            print *, "solve_linear_system failed in compute_boundary_poly (y) at k=", k
            stop
        endif
        
        lag_grid%mk(5:7, k) = xx(1:3)
        
        ! Compute Normal
        if ( k == 1) then
            tp1 =  lag_grid%mk(3, k) ! s=0
            tp2 =  lag_grid%mk(6, k) ! s=0
        else
            tp1 =  2._dp*lag_grid%mk(4, k)*lag_grid%mk(1, k)+lag_grid%mk(3, k)
            tp2 =  2._dp*lag_grid%mk(7, k)*lag_grid%mk(1, k)+lag_grid%mk(6, k)
        endif
        tmp = sqrt(tp1*tp1+tp2*tp2)
        if (tmp < 1.0e-14_dp) then
            print *, 'Error: Zero tangent vector at boundary point k=', k
            print *, 'tp1=', tp1, ' tp2=', tp2
            stop 'Degenerate geometry in compute_boundary_poly'
        endif
        ! Outward normal for Forward (CCW) parametrization: (dy, -dx) -> (tp2, -tp1)
        lag_grid%normal(k, 1) = tp2/tmp
        lag_grid%normal(k, 2) = -tp1/tmp
    enddo

  end subroutine compute_boundary_poly

  subroutine tag_grid(lag_grid, eul_grid, isel)
    type(LagrangianGrid), intent(in) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    integer, intent(in) :: isel
    
    integer :: i,j,nx_local,ny_local,js,kk,side_status
    real(dp) :: x,y,dx,dy,min_dist
    logical :: is_inside
    
    nx_local = eul_grid%nx_grid
    ny_local = eul_grid%ny_grid
    dx = eul_grid%dx
    dy = eul_grid%dy
    
    ! Initialize
    eul_grid%id = 0
    eul_grid%idf = 0
    eul_grid%kaic = -1
    
    do j = 1, ny_local
        do i = 1, nx_local
            ! Optimization: Check if neighbors were tagged in previous step
            ! oid stores previous id (0 or 1/-1?)
            ! In original code: oid stores id (0 or 1)
            ! kdf stores idf (-1 or 1)
            ! If isel != -1, use previous info to skip
            
!            js = abs(eul_grid%oid(i, j)) + abs(eul_grid%oid(i-1, j)) + &
!                 abs(eul_grid%oid(i+1, j)) + abs(eul_grid%oid(i, j-1)) + &
!                 abs(eul_grid%oid(i, j+1))
                 
!            if ((isel /= -1) .and. (js == 0)) then
!                eul_grid%idf(i, j) = eul_grid%chkerc(i, j) ! kdf in original
!                cycle ! goto 678 equivalent
!            endif
            
            x = eul_grid%x_min + (i - xshift) * dx
            y = eul_grid%y_min + (j - yshift) * dy
            
            call classify_nearest_normal_point(x,y,lag_grid%x,lag_grid%y, &
                 lag_grid%normal(:,1),lag_grid%normal(:,2),is_inside,kk, &
                 min_dist,side_status)
            if(side_status/=INTERFACE_SIDE_OK)then
              error stop 'Nearest-normal side classification failed'
            endif
            eul_grid%dmapc(i, j) = sqrt(min_dist)
            eul_grid%kaic(i, j) = kk

            if(is_inside)then
                eul_grid%idf(i, j) = 1
            else
                eul_grid%idf(i, j) = -1
            endif
            
        enddo
    enddo
    
    ! Identify irregular points (interface)
    ! id = 0 if all neighbors have same idf
    ! id = idf if neighbors differ
    do j = 1, ny_local
        do i = 1, nx_local
            js = abs(eul_grid%idf(i, j) - eul_grid%idf(i+1, j)) + &
                 abs(eul_grid%idf(i, j) - eul_grid%idf(i-1, j)) + &
                 abs(eul_grid%idf(i, j) - eul_grid%idf(i, j+1)) + &
                 abs(eul_grid%idf(i, j) - eul_grid%idf(i, j-1))
                 
            if (js == 0) then
                eul_grid%id(i, j) = 0
            else
                eul_grid%id(i, j) = eul_grid%idf(i, j)
            endif
        enddo
    enddo
    
    ! Set kaic to -1 for regular points (id == 0)
    do j = 1, ny_local
        do i = 1, nx_local
            if (eul_grid%id(i, j) == 0) then
                eul_grid%kaic(i, j) = -1
            endif
        enddo
    enddo
    
    ! Boundary correction (first layer 0)
    eul_grid%id(1, :) = 0
    eul_grid%id(nx_local, :) = 0
    eul_grid%id(:, 1) = 0
    eul_grid%id(:, ny_local) = 0
    
    ! Update history arrays
    if (isel == -1) then
        eul_grid%chkero = eul_grid%idf
        eul_grid%chkerc = eul_grid%idf
        eul_grid%oid = eul_grid%id
        eul_grid%qid = eul_grid%id
        
        eul_grid%kaio = eul_grid%kaic
        eul_grid%dmapo = eul_grid%dmapc
    else
        eul_grid%chkero = eul_grid%chkerc
        eul_grid%chkerc = eul_grid%idf
        eul_grid%oid = eul_grid%qid
        eul_grid%qid = eul_grid%id
        
        eul_grid%kaio = eul_grid%kaic
        eul_grid%dmapo = eul_grid%dmapc
    endif
    
    ! Freshly cleared tag
    if (isel /= -1) then
        eul_grid%idn = -eul_grid%chkero + eul_grid%chkerc
    else
        eul_grid%idn = 0
    endif

  end subroutine tag_grid

  subroutine compute_grid_links(lag_grid, eul_grid, isel)
    type(LagrangianGrid), intent(inout) :: lag_grid
    type(EulerianGrid), intent(in) :: eul_grid
    integer, intent(in) :: isel
    
    integer :: k, i, j, im, jm, ik, jk
    real(dp) :: x, y, tx, ty, tdis, tmp, x_c, y_c, x_im, y_jm
    integer :: ico(3,2)
    
    ! Offsets for checking neighbors (0,1), (1,0), (1,1)
    ico(1,1)=0; ico(1,2)=1
    ico(2,1)=1; ico(2,2)=0
    ico(3,1)=1; ico(3,2)=1
    
    do k = 1, lag_grid%npts
        x = lag_grid%x(k)
        y = lag_grid%y(k)
        
        ! 1. Find cell index (i, j) closest to (x, y)
        ! Initial guess: cell containing the point
        i = int((x - eul_grid%x_min) / eul_grid%dx + xshift)
        j = int((y - eul_grid%y_min) / eul_grid%dy + yshift)
        
        ! Calculate cell center coordinates
        x_c = eul_grid%x_min + (real(i, dp) - xshift) * eul_grid%dx
        y_c = eul_grid%y_min + (real(j, dp) - yshift) * eul_grid%dy
        
        tx = x - x_c
        ty = y - y_c
        tdis = tx*tx + ty*ty
        
        jk = 0
        do ik = 1, 3
            im = i + ico(ik, 1)
            jm = j + ico(ik, 2)
            
            x_im = eul_grid%x_min + (real(im, dp) - xshift) * eul_grid%dx
            y_jm = eul_grid%y_min + (real(jm, dp) - yshift) * eul_grid%dy
            
            tmp = (x - x_im)**2 + (y - y_jm)**2
            if (tmp < tdis) then
                tdis = tmp
                jk = ik
            endif
        enddo
        
        if (jk > 0) then
            i = i + ico(jk, 1)
            j = j + ico(jk, 2)
            ! Update tx, ty relative to new closest point
            x_c = eul_grid%x_min + (real(i, dp) - xshift) * eul_grid%dx
            y_c = eul_grid%y_min + (real(j, dp) - yshift) * eul_grid%dy
            tx = x - x_c
            ty = y - y_c
        endif
        
        ! Store center of stencil
        lag_grid%plinkij(k, 1) = i
        lag_grid%plinkij(k, 2) = j
        
        ! 2. Find "corner" point (im, jm)
        ! Initial guess based on quadrant
        im = i + int(sign(1.0_dp, tx))
        jm = j + int(sign(1.0_dp, ty))
        
        ! Check if corner is on the same side of interface
        ! idf stores side (+1 or -1)
        ! We want corner to be on the OPPOSITE side if possible? 
        ! Original code logic:
        ! if (idf(i,j)*idf(im,jm) > 0) then ... try to find opposite side ... endif
        
        if (eul_grid%idf(i, j) * eul_grid%idf(im, jm) > 0) then
            if (abs(tx) > abs(ty)) then
                im = i + int(sign(1.0_dp, tx))
                jm = j - int(sign(1.0_dp, ty))
            else
                im = i - int(sign(1.0_dp, tx))
                jm = j + int(sign(1.0_dp, ty))
            endif
        endif
        
        ! Warning if still on same side (original code prints this)
        if (eul_grid%idf(i, j) * eul_grid%idf(im, jm) > 0) then
             print *, 'corner point is on same side!', isel, k
        endif
        
        lag_grid%plinkij(k, 3) = im
        lag_grid%plinkij(k, 4) = jm
        
    enddo
  end subroutine compute_grid_links

  !-----------------------------------------------------------------------------
  ! Ported Routines for Correction Function
  !-----------------------------------------------------------------------------

  function evalmpoly(s, n, mk, idir, icase, info) result(val)
    integer, intent(in) :: n, idir, icase
    real(dp), intent(in) :: s, mk(7, n)
    integer, intent(out) :: info
    real(dp) :: val
    
    integer :: k
    real(dp) :: ds, s_local, s0, s1, s2
    real(dp) :: c0, c1, c2
    
    info = 0
    ds = tupi / real(n, dp)
    
    ! Find segment k
    ! Assuming mk(1,k) stores s(k) and they are ordered
    ! Simple search
    k = int(s / ds) + 1
    if (k > n) k = n
    if (k < 1) k = 1
    
    ! Adjust k if s is out of [mk(1,k), mk(1,k+1)] due to periodicity or rounding
    ! For now rely on uniform grid mapping
    
    ! mk(1,k) is s_k
    ! mk(2,k), mk(3,k), mk(4,k) are coeffs for x: c0 + c1*s + c2*s^2 ?
    ! No, BdyQuadParametric solves for c0, c1, c2 such that x = c0 + c1*s + c2*s^2
    ! mk(2:4) -> x coeffs
    ! mk(5:7) -> y coeffs
    
    if (idir == 0) then ! x component
        c0 = mk(2, k)
        c1 = mk(3, k)
        c2 = mk(4, k)
    else ! y component
        c0 = mk(5, k)
        c1 = mk(6, k)
        c2 = mk(7, k)
    endif
    
    if (icase == 0) then ! Value
        val = c0 + c1*s + c2*s*s
    elseif (icase == 1) then ! 1st Derivative
        val = c1 + 2.0_dp*c2*s
    elseif (icase == 2) then ! 2nd Derivative
        val = 2.0_dp*c2
    else
        val = 0.0_dp
        info = 1
    endif
    
  end function evalmpoly

  function evalCx(cxcoef, i, x0, center) result(val)
     real(dp), dimension(npts,6), intent(in) :: cxcoef
     integer, intent(in) :: i
     real(dp), dimension(2), intent(in) :: x0, center
     real(dp) :: val
     
     real(dp) :: dx, dy, coe(6)
     
     coe = cxcoef(i, :)
     dx = (x0(1) - center(1))
     dy = (x0(2) - center(2))
     
     ! C(x) = r1 + r2*dx + r3*dy + r4*0.5*dx^2 + r5*0.5*dy^2 + r6*dx*dy
     val = coe(1) + coe(2)*dx + coe(3)*dy + 0.5_dp*dx*dx*coe(4) + &
           0.5_dp*dy*dy*coe(5) + dx*dy*coe(6)
  end function evalCx

  subroutine solve2DCauchy(coef, bdy1, phi_bdy1, bdy2, bdy2_nv, psi_bdy2, bdy_pde, rs_pde, xpos, cxn, dt, time, kappa, vx_pde, vy_pde, steady_state)
      real(dp), dimension(6), intent(out) :: coef
      real(dp), dimension(3,2), intent(in) :: bdy1
      real(dp), dimension(2,2), intent(in) :: bdy2, bdy2_nv
      real(dp), dimension(2), intent(in) :: bdy_pde, psi_bdy2, xpos
      real(dp), dimension(3), intent(in) :: phi_bdy1
      real(dp), intent(in) :: rs_pde, cxn, dt, time, kappa
      real(dp), intent(in) :: vx_pde, vy_pde
      logical, intent(in), optional :: steady_state
      
      real(dp) :: mA(6,6), AF(6,6), rA(6), mR(6), mC(6), xx(6), rcond, &
              ferr(1), berr(1), work(24), dx, dy, vx, vy, nx_val, ny_val, h2
      integer :: iwork(6), ipiv(6), info, i, m
      character(1) :: equed
      
      h2 = hg*hg
      
      ! Setup Matrix mA and RHS rA
      ! Rows 1-3: Value at bdy1 points
      do m = 1, 3
        i = m
        dx = (bdy1(i,1)-xpos(1))/hg
        dy = (bdy1(i,2)-xpos(2))/hg
        mA(i,1) = one
        mA(i,2) = dx
        mA(i,3) = dy
        mA(i,4) = half*dx*dx
        mA(i,5) = half*dy*dy
        mA(i,6) = dx*dy
        rA(i) = phi_bdy1(i) ! Correction matches jump [u] on interface
      enddo
      
      ! Rows 4-5: Normal derivative at bdy2 points
      do m = 1,2
        i = m+3
        dx = (bdy2(m,1) - xpos(1))/hg
        dy = (bdy2(m,2) - xpos(2))/hg
        nx_val = bdy2_nv(m,1)
        ny_val = bdy2_nv(m,2)
        mA(i,1) = zero
        mA(i,2) = nx_val
        mA(i,3) = ny_val
        mA(i,4) = dx*nx_val
        mA(i,5) = dy*ny_val
        mA(i,6) = dy*nx_val+dx*ny_val
        rA(i) = psi_bdy2(m)*hg
      enddo
      
      ! Row 6: PDE at bdy_pde -\Delta u + v.grad u + kappa u = f
      i = 6
      vx = vx_pde
      vy = vy_pde
      dx = (bdy_pde(1)-xpos(1))/hg
      dy = (bdy_pde(2)-xpos(2))/hg
      
! now solve for PDE \partial C/\partial t -Lap u + v.grad u + kappa u = f

      if (.not. present(steady_state) .or. .not. steady_state) then
          mA(i,1) = (kappa+1.0_dp/dt)*h2
          mA(i,2) = vx*hg + (kappa+1.0_dp/dt)*h2*dx
          mA(i,3) = vy*hg + (kappa+1.0_dp/dt)*h2*dy
          mA(i,4) = vx*hg*dx + (kappa+1.0_dp/dt)*0.5_dp*h2*dx*dx - one
          mA(i,5) = vy*hg*dy + (kappa+1.0_dp/dt)*0.5_dp*h2*dy*dy - one
          mA(i,6) = vx*hg*dy + vy*hg*dx + (kappa+1.0_dp/dt)*dx*dy*h2
      else
          mA(i,1) = h2 * kappa ! no h2/dt the time dependent term is on the RHS, this is stationary case
          mA(i,2) = vx*hg + h2*kappa*dx ! stationary
          mA(i,3) = vy*hg + h2*kappa*dy ! stationary
          mA(i,4) = -one + vx*hg*dx + h2*kappa*0.5_dp*dx*dx ! stationary
          mA(i,5) = -one + vy*hg*dy + h2*kappa*0.5_dp*dy*dy ! stationary
          mA(i,6) = vx*hg*dy + vy*hg*dx + h2*kappa*dx*dy
          ! RHS for PDE: F + Cxn/dt for dynamic case
          ! rs_pde passed in is fjmp(i) which is F_jump.
          ! cxn is passed.
      endif
      
      if (present(steady_state) .and. steady_state) then
          rA(i) = rs_pde*h2
      else
          ! Time-dependent case: add C^n/dt term to RHS
          ! cxn is the previous correction evaluated at current point
          rA(i) = rs_pde*h2 + cxn*h2/dt
      endif
      
      ! Solve
      call solve_linear_system(6, mA, rA, xx, rcond, info)
                  
      if (info /= 0) then
          stop
      endif
      
      if (rcond < 1.0e-12_dp) then
          print *, 'Warning: solve2DCauchy matrix ill-conditioned. rcond=', rcond
      endif
      
      coef = xx
      coef(2:3) = coef(2:3)/hg
      coef(4:6) = coef(4:6)/h2
      
  end subroutine solve2DCauchy

  subroutine getCorrection(lag_grid, phi, psi, fjmp, dt, time, Cxcoef, Cxcoefo, &
      kappa,pde_velocity_zero,steady_state,pde_velocity_marker,pde_velocity_scale,status)
    type(LagrangianGrid), intent(in) :: lag_grid
    real(dp), intent(in) :: phi(lag_grid%npts), psi(lag_grid%npts), fjmp(lag_grid%npts)
    real(dp), intent(in) :: dt, time, kappa
    real(dp), intent(out) :: Cxcoef(lag_grid%npts,6) ! (npts, 6)
    real(dp), intent(in), optional :: Cxcoefo(lag_grid%npts,6) ! (npts, 6)
    logical, intent(in), optional :: pde_velocity_zero
    logical, intent(in), optional :: steady_state
    ! Stage 10 may supply the already normalized species velocity at each
    ! marker.  This avoids regenerating or rescaling velocity inside the
    ! correction solve.  Legacy chemical callers omit it and continue through
    ! the accepted frozen-FSI interpolation bridge below.
    real(dp), intent(in), optional :: pde_velocity_marker(:,:)
    real(dp),intent(in),optional::pde_velocity_scale
    integer,intent(out),optional::status
    
    integer :: i, im, ip, info, bridge_status
    real(dp) :: sm, sp, so, sm1, sp1, ds, spar, tmp
    real(dp) :: mA(3,3), pA(3,3), AF(3,3), rA(3), xx(3), mR(3), mC(3), &
            rcond, ferr(1), berr(1), work(12)
    real(dp) :: bdy1_pt(3,2), bdy2_pt(2,2), bdy2_nv(2,2), bdy_pde(2)
    real(dp) :: phi_val(3), psi_val(2), cxn, xt(2), t1, t2
    real(dp) :: vx_pde, vy_pde
    real(dp) :: Cxcoef_trial(lag_grid%npts,6)
    integer :: ipiv(3), iwork(3)
    character(1) :: equed
    logical::explicit_velocity_mode,zero_velocity_mode

    Cxcoef=zero
    Cxcoef_trial=zero
    if(present(status))status=CHEMICAL_FLUX_INVALID
#ifdef SIMCELL_TESTING
    last_correction_marker_velocity=zero
    last_correction_marker_velocity_ready=.false.
#endif
    explicit_velocity_mode=present(pde_velocity_marker)
    zero_velocity_mode=.false.
    if(present(pde_velocity_zero))zero_velocity_mode=pde_velocity_zero
    if(explicit_velocity_mode)then
      if(zero_velocity_mode .or. present(pde_velocity_scale))return
      if(size(pde_velocity_marker,1)/=lag_grid%npts .or. &
          size(pde_velocity_marker,2)/=2)return
      if(.not.all(is_finite_correction_scalar(pde_velocity_marker)))return
    elseif(zero_velocity_mode)then
      if(present(pde_velocity_scale))return
    else
      if(.not.present(pde_velocity_scale))return
      if(.not.is_finite_correction_scalar(pde_velocity_scale))return
      if(pde_velocity_scale<=zero)return
    endif
    
    ds = tupi / real(lag_grid%npts, dp)
    ! Tangential collocation spacing from the moving-domain diffusion scheme.
    spar = 0.5_dp
    
    do i = 1, lag_grid%npts
        so = scc(i)
        
        ! Determine neighbors and s-values
        if (i == 1) then
            im = lag_grid%npts
            ip = 2
            sm = scc(im)-tupi
            sp = scc(ip)
        elseif (i == lag_grid%npts) then
            im = lag_grid%npts - 1
            ip = 1
            sm = scc(im)
            sp = tupi
        else
            im = i - 1
            ip = i + 1
            sm = scc(im)
            sp = scc(ip)
        endif
        
        ! Collocation points for correction value (bdy1)
        sm1 = so - spar*ds
        sp1 = so + spar*ds
        
        bdy1_pt(1,1) = evalmpoly(sm1, lag_grid%npts, lag_grid%mk, 0, 0, info)
        bdy1_pt(1,2) = evalmpoly(sm1, lag_grid%npts, lag_grid%mk, 1, 0, info)
        
        bdy1_pt(2,1) = evalmpoly(so, lag_grid%npts, lag_grid%mk, 0, 0, info)
        bdy1_pt(2,2) = evalmpoly(so, lag_grid%npts, lag_grid%mk, 1, 0, info)
        
        bdy1_pt(3,1) = evalmpoly(sp1, lag_grid%npts, lag_grid%mk, 0, 0, info)
        bdy1_pt(3,2) = evalmpoly(sp1, lag_grid%npts, lag_grid%mk, 1, 0, info)
        
        ! Collocation points for normal derivative (bdy2)
        sm1 = so - spar*ds
        sp1 = so + spar*ds
        
        bdy2_pt(1,1) = evalmpoly(sm1, lag_grid%npts, lag_grid%mk, 0, 0, info)
        bdy2_pt(1,2) = evalmpoly(sm1, lag_grid%npts, lag_grid%mk, 1, 0, info)
        t1 = evalmpoly(sm1, lag_grid%npts, lag_grid%mk, 0, 1, info) ! x'
        t2 = evalmpoly(sm1, lag_grid%npts, lag_grid%mk, 1, 1, info) ! y'
        ! Convert to Normal (y', -x') normalized
        tmp = sqrt(t1**2._dp + t2**2._dp)
        if (tmp < 1.0e-14_dp) then
            print *, 'Warning: Zero tangent at correction point i=', i, ' sm1=', sm1
            tmp = 1.0_dp  ! Fallback to avoid crash
        endif
        bdy2_nv(1,1) = t2 / tmp
        bdy2_nv(1,2) = -t1 / tmp
        
        bdy2_pt(2,1) = evalmpoly(sp1, lag_grid%npts, lag_grid%mk, 0, 0, info)
        bdy2_pt(2,2) = evalmpoly(sp1, lag_grid%npts, lag_grid%mk, 1, 0, info)
        t1 = evalmpoly(sp1, lag_grid%npts, lag_grid%mk, 0, 1, info) ! x'
        t2 = evalmpoly(sp1, lag_grid%npts, lag_grid%mk, 1, 1, info) ! y'
        ! Convert to Normal (y', -x') normalized
        tmp = sqrt(t1**2._dp + t2**2._dp)
        if (tmp < 1.0e-14_dp) then
            print *, 'Warning: Zero tangent at correction point i=', i, ' sp1=', sp1
            tmp = 1.0_dp  ! Fallback to avoid crash
        endif
        bdy2_nv(2,1) = t2 / tmp
        bdy2_nv(2,2) =-t1 / tmp
        
        ! Interpolate psi
        mA(:,1) = one
        mA(1,2) = sm; mA(1,3) = sm*sm
        mA(2,2) = so; mA(2,3) = so*so
        mA(3,2) = sp; mA(3,3) = sp*sp
        
        rA(1) = psi(im); rA(2) = psi(i); rA(3) = psi(ip)
        pA = mA
        
        call solve_linear_system(3, pA, rA, xx, rcond, info)
            
        if (info /= 0) then
            print *, ' getCorrection failed (psi) at i=', i
            stop
        endif
        
        psi_val(1) = xx(1) + xx(2)*sm1 + xx(3)*sm1*sm1
        psi_val(2) = xx(1) + xx(2)*sp1 + xx(3)*sp1*sp1
        
        ! Interpolate phi
        rA(1) = phi(im); rA(2) = phi(i); rA(3) = phi(ip)
        pA = mA
        
        call solve_linear_system(3, pA, rA, xx, rcond, info)
            
        if (info /= 0) then
            print *, ' getCorrection failed (phi) at i=', i
            stop
        endif
        
        phi_val(1) = xx(1) + xx(2)*sm1 + xx(3)*sm1*sm1
        phi_val(2) = xx(1) + xx(2)*so  + xx(3)*so*so
        phi_val(3) = xx(1) + xx(2)*sp1 + xx(3)*sp1*sp1
        
        ! Setup bdy_pde and cxn
        bdy_pde(1) = bdy1_pt(2,1)
        bdy_pde(2) = bdy1_pt(2,2)
        
        cxn = 0.0_dp
        if (present(Cxcoefo)) then
             ! bdy_pde is X^{n+1}, while the stored C^n Taylor
             ! coefficients are centered at the old marker X^n.
             xt(1) = lag_grid%x_old(i)
             xt(2) = lag_grid%y_old(i)
             cxn = evalCx(Cxcoefo, i, bdy_pde, xt)
        endif
        
        if (explicit_velocity_mode) then
            vx_pde = pde_velocity_marker(i,1)
            vy_pde = pde_velocity_marker(i,2)
        elseif (zero_velocity_mode) then
            vx_pde = 0.0_dp
            vy_pde = 0.0_dp
        else
            ! Stage 06: use exactly the accepted FSI velocity frozen before the
            ! chemical solve.  The old implementation regenerated a prescribed
            ! analytic field here, which could disagree with both the FSI state
            ! and the face velocities used by the Cartesian operator.
            call interpolate_fsi_velocity_to_marker(bdy_pde(1), bdy_pde(2), &
                vx_pde, vy_pde, bridge_status)
            if (bridge_status /= VELOCITY_BRIDGE_OK) return
            if(.not.safe_to_scale_correction(vx_pde,pde_velocity_scale) .or. &
                .not.safe_to_scale_correction(vy_pde,pde_velocity_scale))return
            vx_pde=vx_pde*pde_velocity_scale
            vy_pde=vy_pde*pde_velocity_scale
        endif
#ifdef SIMCELL_TESTING
        last_correction_marker_velocity(i,1)=vx_pde
        last_correction_marker_velocity(i,2)=vy_pde
#endif

        call solve2DCauchy(Cxcoef_trial(i,:), bdy1_pt, phi_val, bdy2_pt, bdy2_nv, psi_val, bdy_pde, &
                fjmp(i), bdy_pde, cxn, dt, time, kappa, vx_pde, vy_pde, steady_state)
        
        ! DEBUG: Print info for point 400 (top of circle, y nonzero)
!        if (i == 400) then
!            print '(1x, A)', "  DEBUG getCorrection (i=1):"
!            print '(1x, A, 3F10.4)', "    phi_val(1:3)=", phi_val
!            print '(1x, A, 2F10.4)', "    psi_val(1:2)=", psi_val
!            print '(1x, A, F10.4)', "    fjmp=", fjmp(i)
!            print '(1x, A, 6F10.4)', "    Cxcoef=", Cxcoef(i,:)
!            print '(1x, A, 2F10.4)', "    center (x,y)=", bdy_pde(1), bdy_pde(2)
!        endif
                
    enddo
    Cxcoef=Cxcoef_trial
    if(present(status))status=CHEMICAL_FLUX_OK
#ifdef SIMCELL_TESTING
    last_correction_marker_velocity_ready=.true.
#endif
    return
  end subroutine getCorrection

#ifdef SIMCELL_TESTING
  subroutine get_last_correction_marker_velocity(velocity,status)
    real(dp),intent(out)::velocity(:,:)
    integer,intent(out)::status
    velocity=zero
    status=CHEMICAL_FLUX_INVALID
    if(size(velocity,1)/=npts .or. size(velocity,2)/=2)return
    if(.not.last_correction_marker_velocity_ready)return
    velocity=last_correction_marker_velocity
    status=CHEMICAL_FLUX_OK
  end subroutine get_last_correction_marker_velocity
#endif

  pure elemental logical function is_finite_correction_scalar(value)
    real(dp),intent(in)::value
    integer(int64)::bits
    if(storage_size(value)==64 .and. radix(value)==2 .and. digits(value)==53 .and. &
        minexponent(value)==-1021 .and. maxexponent(value)==1024)then
      bits=transfer(value,bits)
      is_finite_correction_scalar=ibits(bits,52,11)/=int(z'7ff',int64)
    else
      is_finite_correction_scalar=ieee_is_finite(value)
    endif
  end function is_finite_correction_scalar

  pure elemental logical function safe_to_scale_correction(value,scale)
    real(dp),intent(in)::value,scale
    safe_to_scale_correction=.false.
    if(.not.is_finite_correction_scalar(value) .or. &
        .not.is_finite_correction_scalar(scale))return
    if(scale<=zero)return
    if(scale>one)then
      if(abs(value)>huge(one)/scale)return
    endif
    safe_to_scale_correction=.true.
  end function safe_to_scale_correction

  subroutine compute_correction_on_grid(lag_grid, eul_grid, Cxcoef)
    type(LagrangianGrid), intent(in) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    real(dp), intent(in) :: Cxcoef(lag_grid%npts,6)
    
    integer :: i, j, k
    real(dp) :: x, y, xt(2), center(2)
    
    eul_grid%crc = 0.0_dp
    
    ! Loop over Eulerian grid
    do j = 1, eul_grid%ny_grid
        do i = 1, eul_grid%nx_grid
            ! Only compute correction for irregular cells: id = +1 (interior irregular) or id = -1 (exterior irregular)
            if (eul_grid%id(i,j) == 1 .or. eul_grid%id(i,j) == -1) then
                k = eul_grid%kaic(i, j)
                if (k > 0 .and. k <= lag_grid%npts) then
                    xt(1) = eul_grid%x_min + (real(i, dp) - xshift) * eul_grid%dx
                    xt(2) = eul_grid%y_min + (real(j, dp) - yshift) * eul_grid%dy
                    center(1) = lag_grid%x(k)
                    center(2) = lag_grid%y(k)
                    eul_grid%crc(i, j) = evalCx(Cxcoef, k, xt, center)
                else 
                    print *, 'Warning: Invalid kaic at i=', i, ' j=', j, ' kaic=', eul_grid%kaic(i, j)
                    stop
                endif
            endif
        enddo
    enddo
  end subroutine compute_correction_on_grid

  subroutine apply_fresh_cleared_correction(lag_grid, eul_grid, Cxcoefo, field, iside, steady_state, dt, use_old_map, correct_both_sides)
    type(LagrangianGrid), intent(in) :: lag_grid
    type(EulerianGrid), intent(in) :: eul_grid
    real(dp), intent(in) :: Cxcoefo(lag_grid%npts, 6)
    real(dp), intent(inout) :: field(eul_grid%nx_grid, eul_grid%ny_grid)
    integer, intent(in) :: iside
    logical, intent(in) :: steady_state
    real(dp), intent(in) :: dt
    logical, intent(in), optional :: use_old_map
    logical, intent(in), optional :: correct_both_sides
    
    integer :: i, j, k
    real(dp) :: xt(2), center(2), corr
    logical :: both_sides

    ! Backward-compatible default: the accepted chemical path corrects only a
    ! cell newly entering its requested physical side.  Stage 10 actin passes
    ! correct_both_sides=.true. because its full-box auxiliary field represents
    ! both smooth extensions; an interface sweep must therefore add +C^n/dt on
    ! a newly interior cell and -C^n/dt on a newly exterior cell.
    both_sides = .false.
    if (present(correct_both_sides)) both_sides = correct_both_sides
    
    ! Loop over Eulerian grid
    do j = 1, eul_grid%ny_grid
        do i = 1, eul_grid%nx_grid
            ! Check for freshly cleared points
            ! idn = 2 (solid -> fluid) or -2 (fluid -> solid)?
            ! Assuming idn != 0 means status change.
            ! We only care about solid -> fluid (freshly cleared).
            ! Let's assume abs(idn) == 2 implies change.
            ! We need to know which one is "new fluid".
            ! If chkerc (current) is fluid (say 1) and chkero (old) was solid (-1), then idn = -(-1) + 1 = 2.
            ! So idn = 2 is freshly cleared.
            
            !if (abs(eul_grid%idn(i,j)) == 2) then
            if (abs(eul_grid%idn(i,j)) == 2 .and. &
                (both_sides .or. eul_grid%idf(i,j) == iside)) then
                ! The accepted chemical path historically uses kaic.  The
                ! Stage 10 actin transaction passes use_old_map=.true. so an
                ! old Taylor polynomial is paired with its old closest marker.
                if (present(use_old_map)) then
                    if (use_old_map) then
                        k = eul_grid%kaio(i,j)
                    else
                        k = eul_grid%kaic(i,j)
                    endif
                else
                    k = eul_grid%kaic(i,j)
                endif
                if (k > 0 .and. k <= lag_grid%npts) then
                    xt(1) = eul_grid%x_min + (real(i, dp) - xshift) * eul_grid%dx
                    xt(2) = eul_grid%y_min + (real(j, dp) - yshift) * eul_grid%dy
                    center(1) = lag_grid%x_old(k) ! this is for the old correction function in Cxcoefo
                    center(2) = lag_grid%y_old(k) ! this assume old and new pts are not far away from each other
                    
                    ! Compute correction using OLD correction function quadratic polynomial coefficients
                    corr = evalCx(Cxcoefo, k, xt, center)
                    
                    ! Cold is oriented as the iside extension minus its other
                    ! smooth extension.  For the two-sided transaction the
                    ! current tag selects the corresponding +/- history.
                    if (.not. steady_state) then 
                         if (both_sides) then
                             field(i, j) = field(i, j) + &
                                  real(iside*eul_grid%idf(i,j),dp)*corr/dt
                         else
                             field(i, j) = field(i, j) + real(iside,dp)*corr/dt
                         endif
                    endif
                endif
            endif
        enddo
    enddo
  end subroutine apply_fresh_cleared_correction

  subroutine apply_full_correction(lag_grid, eul_grid, Cxcoef, Cxcoefo, iside, steady_state, dt)
    type(LagrangianGrid), intent(in) :: lag_grid
    type(EulerianGrid), intent(inout) :: eul_grid
    real(dp), intent(in) :: Cxcoef(lag_grid%npts, 6)
    real(dp), intent(in), optional :: Cxcoefo(lag_grid%npts, 6)
    integer, intent(in), optional :: iside
    logical, intent(in), optional :: steady_state
    real(dp), intent(in), optional :: dt
    
    integer :: iside_val
    logical :: steady_val
    real(dp) :: dt_val

    if (present(iside)) then
       iside_val = iside
    else
       iside_val = -1 ! Default standard
    endif

    if (present(steady_state)) then
        steady_val = steady_state
    else
        steady_val = .true.
    endif
    
    if (present(dt)) then
        dt_val = dt
    else
        dt_val = 1.0_dp
    endif
    
    ! 1. Compute correction at irregular points (current geometry)
    !    This populates eul_grid%crc
    call compute_correction_on_grid(lag_grid, eul_grid, Cxcoef)
    
    ! NOTE: Freshly cleared correction is NOT applied here.
    ! It should be applied to f_grid (RHS) in compute_H, not to crc.
    
  end subroutine apply_full_correction

  subroutine valatIBpt(avec, ib, cxcoef, lag_grid, eul_grid, cfield, iside)
    real(dp), dimension(6), intent(out) :: avec
    integer, intent(in) :: ib, iside
    type(LagrangianGrid), intent(in) :: lag_grid
    type(EulerianGrid), intent(in) :: eul_grid
    real(dp), dimension(lag_grid%npts, 6), intent(in) :: cxcoef
    real(dp), dimension(-1:eul_grid%nx_grid+1, -1:eul_grid%ny_grid+1), intent(in) :: cfield
    
    real(dp) :: x0, y0, xt, yt
    integer :: i, j, m, cor(6,2)
    real(dp) :: mA(6,6), AF(6,6), rA(6), mR(6), mC(6), xx(6), rcond, &
            ferr(1), berr(1), work(24), dx, dy
    integer :: iwork(6), ipiv(6), info
    character(1) :: equed
    
    real(dp) :: hl, hg2
    real(dp) :: x_center, y_center, tx, ty
    real(dp) :: reconstruction(6),l1_residual,reconstruction_scale
    integer :: i_center, j_center, i_corner, j_corner
    
    x0 = lag_grid%x(ib)
    y0 = lag_grid%y(ib)
    
    hl = eul_grid%dx
    hg2 = hl*hl
    
    ! Use the nearest cell center, its four axial neighbors, and the
    ! diagonal in the marker's quadrant. Selection is independent of side
    ! tags: moving the interface must not replace this full-rank quadratic
    ! stencil with an arbitrary set of six points. Callers keep the interface
    ! inside the supported Cartesian stencil domain.
    i_center = nint((x0 - eul_grid%x_min) / hl + 0.5_dp)
    j_center = nint((y0 - eul_grid%y_min) / hl + 0.5_dp)
    x_center = eul_grid%x_min + (real(i_center, dp) - 0.5_dp) * hl
    y_center = eul_grid%y_min + (real(j_center, dp) - 0.5_dp) * hl
    tx = x0 - x_center
    ty = y0 - y_center

    ! At zero offset, choose the positive direction deterministically.
    i_corner = i_center + 1
    if (tx < 0.0_dp) i_corner = i_center - 1
    j_corner = j_center + 1
    if (ty < 0.0_dp) j_corner = j_center - 1

    ! Preserve the existing row order: up, right, down, left, center, diagonal.
    cor(1,:) = [i_center, j_center + 1]
    cor(2,:) = [i_center + 1, j_center]
    cor(3,:) = [i_center, j_center - 1]
    cor(4,:) = [i_center - 1, j_center]
    cor(5,:) = [i_center, j_center]
    cor(6,:) = [i_corner, j_corner]

    ! =========================================================================
    ! Build interpolation matrix and corrected RHS
    ! =========================================================================
    do m = 1, 6
        i = cor(m,1); j = cor(m,2)
        
        ! Calculate cell center
        xt = eul_grid%x_min + (real(i, dp) - 0.5_dp) * hl
        yt = eul_grid%y_min + (real(j, dp) - 0.5_dp) * hl
        
        dx = xt - x0
        dy = yt - y0
        
        mA(m,1) = one
        mA(m,2) = dx
        mA(m,3) = dy
        mA(m,4) = half*dx*dx
        mA(m,5) = half*dy*dy
        mA(m,6) = dx*dy
        
        ! Side tags determine only the value correction, never the stencil.
        ! The same local jump polynomial applies to every opposite-side point.
        if (eul_grid%idf(i,j) == iside) then
            rA(m) = cfield(i,j)
        else
            rA(m) = cfield(i,j) + real(iside, dp) * evalCx(cxcoef, ib, [xt, yt], [x0, y0])
        endif
    enddo
    
    ! DEBUG: Print stencil info for first few IB points
!    if (ib <= 2) then
!        print *, 'DEBUG valatIBpt ib=', ib, ' iside=', iside, ' x0=', x0, ' y0=', y0
!        do m = 1, 6
!            i = cor(m,1); j = cor(m,2)
!            xt = eul_grid%x_min + (real(i, dp) - 0.5_dp) * hl
!            yt = eul_grid%y_min + (real(j, dp) - 0.5_dp) * hl
!            print *, '  Point ', m, ': (', i, ',', j, ') at (', xt, ',', yt, ')'
!            print *, '    idf=', eul_grid%idf(i,j), ' id=', eul_grid%id(i,j), ' field=', cfield(i,j), ' rA=', rA(m)
!        enddo
!    endif
    
    ! Solve system
    call solve_linear_system(6, mA, rA, xx, rcond, info)
    fixed_stencil_calls=fixed_stencil_calls+1_int64
    if(info/=0)then
      fixed_stencil_rank_failures=fixed_stencil_rank_failures+1_int64
      fixed_stencil_reconstruction_failures= &
           fixed_stencil_reconstruction_failures+1_int64
    else
      if(ieee_is_finite(rcond) .and. rcond>zero) &
           fixed_stencil_min_rcond=min(fixed_stencil_min_rcond,rcond)
      if(.not.ieee_is_finite(rcond) .or. &
           rcond<=100.0_dp*epsilon(one)) &
           fixed_stencil_near_rank_deficient= &
           fixed_stencil_near_rank_deficient+1_int64
      reconstruction=matmul(mA,xx)-rA
      l1_residual=sum(abs(reconstruction))
      reconstruction_scale=max(one,sum(abs(rA)), &
           sum(matmul(abs(mA),abs(xx))))
      if(ieee_is_finite(l1_residual) .and. &
           ieee_is_finite(reconstruction_scale))then
        fixed_stencil_max_l1_residual=max( &
             fixed_stencil_max_l1_residual,l1_residual)
        fixed_stencil_max_scaled_l1_residual=max( &
             fixed_stencil_max_scaled_l1_residual, &
             l1_residual/reconstruction_scale)
        if(l1_residual>4096.0_dp*epsilon(one)*reconstruction_scale) &
             fixed_stencil_reconstruction_failures= &
             fixed_stencil_reconstruction_failures+1_int64
      else
        fixed_stencil_reconstruction_failures= &
             fixed_stencil_reconstruction_failures+1_int64
      endif
    endif
            
    if (info /= 0) then
        ! Matrix is singular or ill-conditioned (duplicate points near boundary)
        ! Fallback: return simple interpolation at IB point
        print *, 'WARNING: valatIBpt solve_linear_system failed (info=', info, ') ib=', ib, &
                ' - using field value at closest point'
        
        ! Use field value at center point (cor(5,:)) as fallback
        i = cor(5,1); j = cor(5,2)
        xt = eul_grid%x_min + (real(i, dp) - 0.5_dp) * hl
        yt = eul_grid%y_min + (real(j, dp) - 0.5_dp) * hl
        
        ! Set output: value at center, zero gradient
        avec(1) = cfield(i,j)
        avec(2) = 0.0_dp  ! du/dx
        avec(3) = 0.0_dp  ! du/dy
        avec(4) = 0.0_dp  ! d2u/dx2
        avec(5) = 0.0_dp  ! d2u/dy2
        avec(6) = 0.0_dp  ! d2u/dxdy
        return
        stop
    endif
    
    avec = xx
    
  end subroutine valatIBpt

  !============================================================================
  ! interp_field_to_lagpoint: Safe interpolation from Eulerian field to Lagrangian point
  !
  ! Purpose: Interpolate a field (from previous time level) to a Lagrangian point
  !          at the NEW position, using Eulerian grid points that are on the
  !          same side in the CURRENT geometry.
  !
  ! This is needed for coupled problems where:
  ! - Interface has moved from old to new position
  ! - We need BC values at new interface positions using old field values
  !
  ! Algorithm:
  ! 1. Use CURRENT geometry (idf) to select points on side iside (4x4 patch)
  ! 2. Check if each point is "freshly cleared" (chkero != idf)
  ! 3. If freshly cleared, apply correction using PRE-COMPUTED crc:
  !    - For interior (iside=1):  u_in = u_out + crc
  !    - For exterior (iside=-1): u_out = u_in - crc
  ! 4. Build least-squares polynomial fit using (corrected) values
  !============================================================================
  subroutine interp_field_to_lagpoint(x0, y0, eul_grid, cfield, iside, val_out, grad_out)
    real(dp), intent(in) :: x0, y0                                    ! Target Lagrangian point
    type(EulerianGrid), intent(in) :: eul_grid                        ! Eulerian grid (with chkero and crc from previous step)
    real(dp), dimension(-1:,:), intent(in) :: cfield                  ! Field to interpolate
    integer, intent(in) :: iside                                      ! Which side: +1 = interior, -1 = exterior
    real(dp), intent(out) :: val_out                                  ! Interpolated value
    real(dp), dimension(2), intent(out), optional :: grad_out         ! Interpolated gradient (du/dx, du/dy)
    
    ! Local variables
    integer :: ii, jj, i_center, j_center, m
    integer :: i_lo, i_hi, j_lo, j_hi
    real(dp) :: hl, dx, dy, corr_val
    logical :: is_freshly_cleared
    
    ! Storage for candidate points (up to 16 points in 4x4 stencil)
    integer, parameter :: MAX_CANDIDATES = 16
    integer :: cand_i(MAX_CANDIDATES), cand_j(MAX_CANDIDATES)
    real(dp) :: cand_x(MAX_CANDIDATES), cand_y(MAX_CANDIDATES)
    real(dp) :: cand_val(MAX_CANDIDATES)
    integer :: n_cand
    
    ! Least squares matrix (6 coefficients for quadratic: 1, x, y, x^2, y^2, xy)
    real(dp) :: mA(MAX_CANDIDATES, 6), rhs_vec(MAX_CANDIDATES)
    real(dp) :: ATA(6,6), ATb(6), coef(6), rcond
    integer :: info
    
    hl = eul_grid%dx
    
    ! =========================================================================
    ! STEP 1: Find cell center closest to the Lagrangian point
    ! =========================================================================
    i_center = nint((x0 - eul_grid%x_min) / hl + 0.5_dp)
    j_center = nint((y0 - eul_grid%y_min) / hl + 0.5_dp)
    
    ! Clamp to valid range (need room for 4x4 stencil)
    i_center = max(2, min(i_center, eul_grid%nx_grid - 1))
    j_center = max(2, min(j_center, eul_grid%ny_grid - 1))
    
    ! =========================================================================
    ! STEP 2: Collect candidate points from 4x4 stencil
    !         Select points where idf == iside (same side in CURRENT geometry)
    !         Apply correction to freshly cleared points using stored crc
    ! =========================================================================
    ! 4x4 patch: from (i_center-1, j_center-1) to (i_center+2, j_center+2)
    i_lo = i_center - 1
    i_hi = i_center + 2
    j_lo = j_center - 1
    j_hi = j_center + 2
    
    ! Clamp to valid range
    i_lo = max(1, i_lo); i_hi = min(eul_grid%nx_grid, i_hi)
    j_lo = max(1, j_lo); j_hi = min(eul_grid%ny_grid, j_hi)
    
    n_cand = 0
    do jj = j_lo, j_hi
        do ii = i_lo, i_hi
            ! Check if this point is on the correct side in CURRENT geometry
            if (eul_grid%idf(ii, jj) == iside) then
                n_cand = n_cand + 1
                if (n_cand > MAX_CANDIDATES) exit
                
                cand_i(n_cand) = ii
                cand_j(n_cand) = jj
                cand_x(n_cand) = eul_grid%x_min + (real(ii, dp) - 0.5_dp) * hl
                cand_y(n_cand) = eul_grid%y_min + (real(jj, dp) - 0.5_dp) * hl
                
                ! Check if this is a freshly cleared point
                ! Freshly cleared = was on OTHER side at previous time level
                is_freshly_cleared = (eul_grid%chkero(ii, jj) /= iside)
                
                if (is_freshly_cleared) then
                    ! This point was on the other side at previous time
                    ! The stored field value is from the other side's solution
                    ! Use the PRE-COMPUTED crc value from the previous step
                    ! crc = [u] = u_in - u_out was computed at this grid point last step
                    corr_val = eul_grid%crc(ii, jj)
                    
                    ! Apply correction: 
                    ! For interior (iside=1):  u_in = u_out + crc (add correction)
                    ! For exterior (iside=-1): u_out = u_in - crc (subtract correction)
                    cand_val(n_cand) = cfield(ii, jj) + real(iside, dp) * corr_val
                else
                    ! Point was on our side at previous time - use value directly
                    cand_val(n_cand) = cfield(ii, jj)
                endif
            endif
        enddo
        if (n_cand > MAX_CANDIDATES) exit
    enddo
    
    ! =========================================================================
    ! STEP 3: Check if we have enough points for interpolation
    ! =========================================================================
    if (n_cand < 6) then
        ! Fallback: expand search to 6x6 stencil
        i_lo = max(1, i_center - 2)
        i_hi = min(eul_grid%nx_grid, i_center + 3)
        j_lo = max(1, j_center - 2)
        j_hi = min(eul_grid%ny_grid, j_center + 3)
        
        n_cand = 0
        do jj = j_lo, j_hi
            do ii = i_lo, i_hi
                if (eul_grid%idf(ii, jj) == iside) then
                    n_cand = n_cand + 1
                    if (n_cand > MAX_CANDIDATES) exit
                    
                    cand_i(n_cand) = ii
                    cand_j(n_cand) = jj
                    cand_x(n_cand) = eul_grid%x_min + (real(ii, dp) - 0.5_dp) * hl
                    cand_y(n_cand) = eul_grid%y_min + (real(jj, dp) - 0.5_dp) * hl
                    
                    is_freshly_cleared = (eul_grid%chkero(ii, jj) /= iside)
                    
                    if (is_freshly_cleared) then
                        corr_val = eul_grid%crc(ii, jj)
                        cand_val(n_cand) = cfield(ii, jj) + real(iside, dp) * corr_val
                    else
                        cand_val(n_cand) = cfield(ii, jj)
                    endif
                endif
            enddo
            if (n_cand > MAX_CANDIDATES) exit
        enddo
    endif
    
    if (n_cand < 3) then
        ! Emergency fallback: use nearest point value
        print *, 'WARNING: interp_field_to_lagpoint - only', n_cand, &
                 'points found near (', x0, ',', y0, ')'
        if (n_cand >= 1) then
            val_out = cand_val(1)
        else
            val_out = 0.0_dp
        endif
        if (present(grad_out)) grad_out = 0.0_dp
        return
    endif
    
    ! =========================================================================
    ! STEP 4: Build least-squares system for polynomial fit
    !         f(x,y) = c1 + c2*(x-x0) + c3*(y-y0) + c4*(x-x0)^2/2 + c5*(y-y0)^2/2 + c6*(x-x0)*(y-y0)
    ! =========================================================================
    mA = 0.0_dp
    rhs_vec = 0.0_dp
    
    do m = 1, n_cand
        dx = cand_x(m) - x0
        dy = cand_y(m) - y0
        
        mA(m, 1) = 1.0_dp
        mA(m, 2) = dx
        mA(m, 3) = dy
        mA(m, 4) = 0.5_dp * dx * dx
        mA(m, 5) = 0.5_dp * dy * dy
        mA(m, 6) = dx * dy
        
        rhs_vec(m) = cand_val(m)
    enddo
    
    ! =========================================================================
    ! STEP 5: Solve least-squares system A^T A c = A^T b
    ! =========================================================================
    ! Form A^T A (6x6)
    ATA = 0.0_dp
    ATb = 0.0_dp
    
    do m = 1, n_cand
        do ii = 1, 6
            ATb(ii) = ATb(ii) + mA(m, ii) * rhs_vec(m)
            do jj = 1, 6
                ATA(ii, jj) = ATA(ii, jj) + mA(m, ii) * mA(m, jj)
            enddo
        enddo
    enddo
    
    ! Solve the normal equations
    call solve_linear_system(6, ATA, ATb, coef, rcond, info)
    
    if (info /= 0) then
        ! Fallback: use simple average
        print *, 'WARNING: interp_field_to_lagpoint - least squares failed, using average'
        val_out = sum(cand_val(1:n_cand)) / real(n_cand, dp)
        if (present(grad_out)) grad_out = 0.0_dp
        return
    endif
    
    ! =========================================================================
    ! STEP 6: Evaluate polynomial at x0, y0 (which gives c1 directly)
    ! =========================================================================
    val_out = coef(1)  ! At (x0, y0), only the constant term survives
    
    if (present(grad_out)) then
        grad_out(1) = coef(2)  ! du/dx at (x0, y0)
        grad_out(2) = coef(3)  ! du/dy at (x0, y0)
    endif
    
  end subroutine interp_field_to_lagpoint

  !---------------------------------------------------------------------------
  ! Simple bilinear interpolation fallback
  !---------------------------------------------------------------------------
  subroutine interp_field_bilinear(x0, y0, eul_grid, cfield, iside, val_out)
    real(dp), intent(in) :: x0, y0
    type(EulerianGrid), intent(in) :: eul_grid
    real(dp), dimension(-1:,:), intent(in) :: cfield
    integer, intent(in) :: iside
    real(dp), intent(out) :: val_out
    
    integer :: i, j, ii, jj
    real(dp) :: hl, u, v, x_center, y_center
    real(dp) :: w_sum, val_sum, w
    
    hl = eul_grid%dx
    
    ! Find closest cell center
    i = nint((x0 - eul_grid%x_min) / hl + 0.5_dp)
    j = nint((y0 - eul_grid%y_min) / hl + 0.5_dp)
    i = max(2, min(i, eul_grid%nx_grid - 1))
    j = max(2, min(j, eul_grid%ny_grid - 1))
    
    ! Weighted average of nearby points on correct side
    w_sum = 0.0_dp
    val_sum = 0.0_dp
    do jj = j-1, j+1
        do ii = i-1, i+1
            if (ii < 1 .or. ii > eul_grid%nx_grid) cycle
            if (jj < 1 .or. jj > eul_grid%ny_grid) cycle
            if (eul_grid%idf(ii, jj) /= iside) cycle
            
            x_center = eul_grid%x_min + (real(ii, dp) - 0.5_dp) * hl
            y_center = eul_grid%y_min + (real(jj, dp) - 0.5_dp) * hl
            w = 1.0_dp / (((x0 - x_center)**2 + (y0 - y_center)**2) + 1.0e-10_dp)
            
            ! Apply freshly-cleared correction if needed
            if (abs(eul_grid%idn(ii, jj)) == 2) then
                val_sum = val_sum + w * (cfield(ii, jj) - real(iside, dp) * eul_grid%crc(ii, jj))
            else
                val_sum = val_sum + w * cfield(ii, jj)
            endif
            w_sum = w_sum + w
        enddo
    enddo
    
    if (w_sum > 0.0_dp) then
        val_out = val_sum / w_sum
    else
        ! Last resort: use center value
        val_out = cfield(i, j)
    endif
    
  end subroutine interp_field_bilinear

end module geometry_mod
