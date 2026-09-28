/*=============================================================================
*
*   Filename : AdvecDiffMG.cpp
*   Creator : Han Zhou
*   Date : 11/16/25
*   Description : This code solves the following equation:
*       -Delta u + div(v u) + kappa u = f,
*       in the domain (low[0], high[0]) x (low[1], high[1])
*       with periodic BC in x and homogeneous total no-flux BC on the
*       horizontal (south/north) boundaries.
*       We use a uniform Cartesian grid with (Nx x Ny) cells with
*       cell-centered discretization.
*       The velocity field is given at edge centers.
*
=============================================================================*/

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "Variables.h"
#include "MultiLevelSolver.h"

void vec2mat(const double *v, int len, MatrixXd &mat, int row, int col)
{
  assert(len == row*col);
  for(int j = 0; j < col; j++){
    for(int i = 0; i < row; i++){
      mat[i][j] = v[i + row*j];
    }
  }
}
void mat2vec(const MatrixXd &mat, int row, int col, double *v, int len)
{
  assert(len == row*col);
  for(int j = 0; j < col; j++){
    for(int i = 0; i < row; i++){
      v[i + row*j] = mat[i][j];
    }
  }
}

// Solve with fixed periodic-x / homogeneous-total-no-flux-y semantics.
extern "C"
int solveAdvecDiffEqn(double Lx, double Ly,
                      double Hx, double Hy,
                      int Nx, int Ny,
                      const double *vx,   // (Nx+1)*Ny, east/west edge ctr
                      const double *vy,   // Nx*(Ny+1), north/south edge ctr
                      double kappa,       // reaction coefficient
                      const double *f,    // Nx*Ny, RHS at cell center
                      double rel_tol,     // e.g. 1e-8
                      double abs_tol,     // e.g. 1e-10
                      int maxitrn,        // max iteration number, e.g. 20
                      double *u)          // Nx*Ny, solution at cell center
{
  double low[2] = {Lx, Ly};
  double high[2] = {Hx, Hy};

  MatrixXd vx_mat(Nx+1, Ny, 0.0),
           vy_mat(Nx, Ny+1, 0.0),
           f_mat(Nx, Ny, 0.0),
           u_mat(Nx, Ny, 0.0);

  vec2mat(vx, (Nx+1)*Ny, vx_mat, Nx+1, Ny);
  vec2mat(vy, Nx*(Ny+1), vy_mat, Nx, Ny+1);
  vec2mat(f, Nx*Ny, f_mat, Nx, Ny);
  vec2mat(u, Nx*Ny, u_mat, Nx, Ny);

  MultiLevelSolver solver(low, high, Nx, Ny, vx_mat, vy_mat, kappa);

  int itrn = solver.solveWithMultigrid(f_mat, u_mat, rel_tol, abs_tol, maxitrn);

  mat2vec(u_mat, Nx, Ny, u, Nx*Ny);
  return itrn;
}

// Solve the interface problem with periodic x and homogeneous total no-flux
// at the south/north walls.  The north/south arrays remain in the unchanged
// C ABI but are unused because this implementation fixes those wall data to
// zero, matching the original solver's behavior.
extern "C"
int solveAdvecDiffInterfaceProblem(double Lx, double Ly,
                                   double Hx, double Hy,
                                   int Nx, int Ny,
                                   const double *vx,
                                   const double *vy,
                                   double kappa,
                                   const double *f,
                                   const double *bc_n,
                                   const double *bc_s,
                                   const double *mask,
                                   const double *crc,
                                   double rel_tol,
                                   double abs_tol,
                                   int maxitrn,
                                   double *u)
{
  (void)bc_n;
  (void)bc_s;

  double low[2] = {Lx, Ly};
  double high[2] = {Hx, Hy};

  MatrixXd vx_mat(Nx+1, Ny, 0.0),
           vy_mat(Nx, Ny+1, 0.0),
           f_mat(Nx, Ny, 0.0),
           u_mat(Nx, Ny, 0.0),
           mask_mat(Nx, Ny, 0.0),
           crc_mat(Nx, Ny, 0.0);

  vec2mat(vx, (Nx+1)*Ny, vx_mat, Nx+1, Ny);
  vec2mat(vy, Nx*(Ny+1), vy_mat, Nx, Ny+1);
  vec2mat(f, Nx*Ny, f_mat, Nx, Ny);
  vec2mat(u, Nx*Ny, u_mat, Nx, Ny);
  vec2mat(mask, Nx*Ny, mask_mat, Nx, Ny);
  vec2mat(crc, Nx*Ny, crc_mat, Nx, Ny);

  MultiLevelSolver solver(low, high, Nx, Ny, vx_mat, vy_mat, kappa);

  solver.correctRightHandSide(mask_mat, crc_mat, f_mat);

  int itrn = solver.solveWithMultigrid(f_mat, u_mat, rel_tol, abs_tol, maxitrn);

  mat2vec(u_mat, Nx, Ny, u, Nx*Ny);
  return itrn;
}
