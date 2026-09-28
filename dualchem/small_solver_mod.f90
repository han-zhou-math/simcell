module small_solver_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: solve_linear_system, qr_decompose, solve_qr_system
  public :: report_small_solver_audit
  
  integer, parameter :: dp = selected_real_kind(p = 15, r = 300)
  integer(kind=8),save :: audit_solve_count=0_8
  integer(kind=8),save :: audit_failure_count=0_8
  integer(kind=8),save :: audit_nonfinite_rcond_count=0_8
  integer(kind=8),save :: audit_small_rcond_count=0_8
  real(dp),save :: audit_min_positive_rcond=huge(1.0_dp)

contains

  ! Wrapper that now uses QR by default (or we can switch)
  subroutine solve_linear_system(n, A_in, b_in, x, rcond, info)
    integer, intent(in) :: n
    real(dp), intent(in) :: A_in(n,n)
    real(dp), intent(in) :: b_in(n)
    real(dp), intent(out) :: x(n)
    real(dp), intent(out) :: rcond
    integer, intent(out) :: info
    
    real(dp) :: Q(n,n), R(n,n)
    
    ! Call QR solver
    call qr_decompose(n, A_in, Q, R, rcond, info)
    audit_solve_count=audit_solve_count+1_8
    if(info/=0)then
      audit_failure_count=audit_failure_count+1_8
    elseif(.not.ieee_is_finite(rcond))then
      audit_nonfinite_rcond_count=audit_nonfinite_rcond_count+1_8
    else
      if(rcond>0.0_dp)audit_min_positive_rcond=min( &
           audit_min_positive_rcond,rcond)
      if(rcond<=100.0_dp*epsilon(1.0_dp)) &
           audit_small_rcond_count=audit_small_rcond_count+1_8
    endif
    
    if (info /= 0) then
         return
    endif
    
    call solve_qr_system(n, Q, R, b_in, x)

  end subroutine solve_linear_system

  subroutine report_small_solver_audit()
    real(dp)::reported_minimum

    reported_minimum=audit_min_positive_rcond
    if(audit_solve_count==0_8 .or. &
         reported_minimum==huge(1.0_dp))reported_minimum=0.0_dp
    write(*,'(a,1x,4(i0,1x),es24.16)') 'FIG2_STENCIL_SOLVER_AUDIT', &
         audit_solve_count,audit_failure_count, &
         audit_nonfinite_rcond_count,audit_small_rcond_count, &
         reported_minimum
  end subroutine report_small_solver_audit

  ! Householder QR Decomposition
  ! A = Q * R
  ! Input: A (n,n)
  ! Output: Q (n,n) - Orthogonal matrix
  !         R (n,n) - Upper triangular matrix
  !         rcond   - Estimate of condition number (min_diag / max_diag of R)
  !         info    - 0 if successful, k if R(k,k) is singular
  subroutine qr_decompose(n, A_in, Q, R, rcond, info)
    integer, intent(in) :: n
    real(dp), intent(in) :: A_in(n,n)
    real(dp), intent(out) :: Q(n,n), R(n,n)
    real(dp), intent(out) :: rcond
    integer, intent(out) :: info
    
    integer :: k, i, j
    real(dp) :: v(n), norm_x, alpha, u_norm_sq, tau, dot_prod
    real(dp) :: H(n,n), temp_mat(n,n)
    real(dp) :: min_diag, max_diag
    
    R = A_in
    Q = 0.0_dp
    do i = 1, n
        Q(i,i) = 1.0_dp
    enddo
    
    info = 0
    rcond = 0.0_dp
    
    do k = 1, n-1
        ! Compute Householder vector v for column k below diagonal
        norm_x = 0.0_dp
        do i = k, n
            norm_x = norm_x + R(i,k)**2
        enddo
        norm_x = sqrt(norm_x)
        
        if (norm_x < 1.0e-30_dp) then
            ! Column is already zero? Continue
            cycle
        endif
        
        alpha = -sign(norm_x, R(k,k))
        
        v = 0.0_dp
        v(k) = R(k,k) - alpha
        do i = k+1, n
            v(i) = R(i,k)
        enddo
        
        u_norm_sq = 0.0_dp
        do i = k, n
            u_norm_sq = u_norm_sq + v(i)**2
        enddo
        
        if (u_norm_sq < 1.0e-30_dp) cycle
        
        ! tau = 2 / (v'v)
        tau = 2.0_dp / u_norm_sq
        
        ! Apply H = I - tau * v * v' to R from left: R = H * R
        ! R(k:n, k:n) = R(k:n, k:n) - tau * v * (v' * R(k:n, k:n))
        
        ! w = v' * R
        ! w(j) = sum_i v(i) * R(i,j)
        do j = k, n
            dot_prod = 0.0_dp
            do i = k, n
                dot_prod = dot_prod + v(i) * R(i,j)
            enddo
            
            do i = k, n
                R(i,j) = R(i,j) - tau * v(i) * dot_prod
            enddo
        enddo
        
        ! Apply H to Q from right: Q = Q * H' = Q * H
        ! Q(1:n, k:n) = Q(1:n, k:n) - tau * (Q(1:n, k:n) * v) * v'
        
        ! w = Q * v
        do i = 1, n
            dot_prod = 0.0_dp
            do j = k, n
                dot_prod = dot_prod + Q(i,j) * v(j)
            enddo
            
            do j = k, n
                Q(i,j) = Q(i,j) - tau * dot_prod * v(j)
            enddo
        enddo
    enddo
    
    ! Check singularity
    max_diag = 0.0_dp
    min_diag = huge(1.0_dp)
    
    do k = 1, n
        if (abs(R(k,k)) < 1.0e-30_dp) then
            info = k
            return
        endif
        max_diag = max(max_diag, abs(R(k,k)))
        min_diag = min(min_diag, abs(R(k,k)))
    enddo
    
    if (max_diag > 0.0_dp) then
        rcond = min_diag / max_diag
    endif

  end subroutine qr_decompose

  ! Solve R x = Q' b
  subroutine solve_qr_system(n, Q, R, b_in, x)
    integer, intent(in) :: n
    real(dp), intent(in) :: Q(n,n), R(n,n), b_in(n)
    real(dp), intent(out) :: x(n)
    
    real(dp) :: y(n), dot_prod
    integer :: i, j
    
    ! y = Q' * b
    do i = 1, n
        dot_prod = 0.0_dp
        do j = 1, n
            dot_prod = dot_prod + Q(j,i) * b_in(j)
        enddo
        y(i) = dot_prod
    enddo
    
    ! Solve R x = y (Back substitution)
    x(n) = y(n) / R(n,n)
    do i = n-1, 1, -1
        dot_prod = y(i)
        do j = i+1, n
            dot_prod = dot_prod - R(i,j) * x(j)
        enddo
        x(i) = dot_prod / R(i,i)
    enddo
    
  end subroutine solve_qr_system

end module small_solver_mod
