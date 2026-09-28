module advec_diff_solver_mod
  use iso_c_binding
  use parameters, only: dp, dualchem_mg_rtol, dualchem_mg_atol, &
      dualchem_mg_max_iterations
  implicit none
  private
  public :: solve_advec_diff

  interface
    function c_solve_advec_diff_interface(Lx, Ly, Hx, Hy, Nx, Ny, &
                                vx, vy, kappa, f, &
                                bc_n, bc_s, &
                                mask, crc, &
                                rel_tol, abs_tol, maxitrn, u) bind(C, name="solveAdvecDiffInterfaceProblem")
      import :: c_double, c_int
      real(c_double), value :: Lx, Ly, Hx, Hy
      integer(c_int), value :: Nx, Ny
      real(c_double), intent(in) :: vx(*), vy(*)
      real(c_double), value :: kappa
      real(c_double), intent(in) :: f(*)
      real(c_double), intent(in) :: bc_n(*), bc_s(*)
      real(c_double), intent(in) :: mask(*), crc(*)
      real(c_double), value :: rel_tol, abs_tol
      integer(c_int), value :: maxitrn
      real(c_double), intent(inout) :: u(*)
      integer(c_int) :: c_solve_advec_diff_interface
    end function c_solve_advec_diff_interface
  end interface

contains

  ! Wrapper routine to provide a cleaner Fortran API and handle any potential type/shape mismatches
  ! Note: Fortran arrays are column-major. The C++ solver's vec2mat function expects
  ! input that corresponds to column-major flattening (v[i + row*j]), so direct passing works.
  ! Pass the first storage element of each array so the C binding receives a flat base address.
  subroutine solve_advec_diff(Lx, Ly, Hx, Hy, Nx, Ny, vx, vy, kappa, f, bc_n, bc_s, mask, crc, u, iter, converged)
    real(dp), intent(in) :: Lx, Ly, Hx, Hy
    integer, intent(in) :: Nx, Ny
    real(dp), intent(in) :: vx(Nx+1, Ny), vy(Nx, Ny+1)
    real(dp), intent(in) :: kappa
    real(dp), intent(in) :: f(Nx, Ny)
    real(dp), intent(in) :: bc_n(Nx), bc_s(Nx)
    real(dp), intent(in) :: mask(Nx, Ny), crc(Nx, Ny)
    real(dp), intent(inout) :: u(Nx, Ny)
    integer, intent(out) :: iter
    logical, intent(out) :: converged
    
    integer(c_int) :: c_iter
    
    ! Call the C function using the base storage address for each Fortran array.
    c_iter = c_solve_advec_diff_interface( &
        real(Lx, c_double), real(Ly, c_double), &
        real(Hx, c_double), real(Hy, c_double), &
        int(Nx, c_int), int(Ny, c_int), &
      vx(1,1), vy(1,1), &
        real(kappa, c_double), &
      f(1,1), &
      bc_n(1), bc_s(1), &
      mask(1,1), crc(1,1), &
        real(dualchem_mg_rtol,c_double), real(dualchem_mg_atol,c_double), &
        int(dualchem_mg_max_iterations,c_int), &
      u(1,1))
        
    iter = int(c_iter)
    converged = iter /= dualchem_mg_max_iterations
    
  end subroutine solve_advec_diff

end module advec_diff_solver_mod
