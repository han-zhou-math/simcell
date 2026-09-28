/*=============================================================================
*
*   Filename : MultiLevelSolver.h
*   Creator : Han Zhou
*   Date : 11/14/25
*   Description :
*
=============================================================================*/

#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "Variables.h"

class MultiLevelSolver
{
  private:

    double _low[2]  = {0.0,0.0};
    double _high[2] = {1.0,1.0};

    MatrixXd _vx;
    MatrixXd _vy;

    double _kappa = 0.0;

    int _Nx = 0;
    int _Ny = 0;
    double _dx = 0.0;
    double _dy = 0.0;
    double _dx2 = 0.0;
    double _dy2 = 0.0;
    double _one_dx = 0.0;
    double _one_dy = 0.0;
    double _one_dx2 = 0.0;
    double _one_dy2 = 0.0;

    Matrix<Array<double,5>> _fds_coeff;
    MatrixXd _one_diag;

    MultiLevelSolver *finer = 0;
    MultiLevelSolver *coarser = 0;

  public:

    MultiLevelSolver(double low[2], double high[2],
                     int Nx, int Ny,
                     const MatrixXd &vx,
                     const MatrixXd &vy,
                     double kappa)
    :_Nx(Nx), _Ny(Ny), _kappa(kappa)
    {
      assert(isPowerTwo(Nx));
      assert(isPowerTwo(Ny));

      for(int i = 0; i < 2; i++){
        _low[i] = low[i];
        _high[i] = high[i];
      }

      _dx = (high[0]-low[0])/Nx;
      _dy = (high[1]-low[1])/Ny;
      _dx2 = _dx*_dx;
      _dy2 = _dy*_dy;
      _one_dx = 1.0/_dx;
      _one_dy = 1.0/_dy;
      _one_dx2 = 1.0/_dx2;
      _one_dy2 = 1.0/_dy2;

      _vx = vx;
      _vy = vy;

      _one_diag.reallocate(Nx, Ny);
      _one_diag.fill(0.0);
      _fds_coeff.reallocate(Nx, Ny);
      _fds_coeff.fill(0.0);
      getFDSCoeff(_fds_coeff);
      getDiagInv(_one_diag);

      if (_Nx >= 2 && _Ny >= 2) {

        int Nx_2 = Nx/2;
        int Ny_2 = Ny/2;

        MatrixXd vx_2(Nx_2+1, Ny_2, 0.0), vy_2(Nx_2, Ny_2+1, 0.0);

        for(int i = 0, i1 = 0; i <= Nx_2; i++, i1 += 2){
          for(int j = 0, j1 = 0; j < Ny_2; j++, j1 += 2){
            vx_2[i][j] = 0.5*(vx[i1][j1] + vx[i1][j1+1]);
          }
        }
        // The final x face is the periodic duplicate of the first face.
        // Assign it explicitly at every MG level so restriction cannot
        // introduce a seam mismatch through independent roundoff.
        for(int j = 0; j < Ny_2; j++){
          vx_2[Nx_2][j] = vx_2[0][j];
        }
        for(int i = 0, i1 = 0; i < Nx_2; i++, i1 += 2){
          for(int j = 0, j1 = 0; j <= Ny_2; j++, j1 += 2){
            vy_2[i][j] = 0.5*(vy[i1][j1] + vy[i1+1][j1]);
          }
        }

        coarser = new MultiLevelSolver(low, high, Nx_2, Ny_2, vx_2, vy_2, kappa);
        coarser->finer = this;
      }
    }

    virtual ~MultiLevelSolver(void)
    {
      finer = 0;
      if (coarser != 0) {
        delete coarser;
        coarser = 0;
      }
    }

    bool isPowerTwo(int k) const
    {
      int p = static_cast<int>(log(k)/log(2.0) + 0.5);
      return k == 1<<p;
    }

    // All stencil consumers share these accessors.  At the 1-cell coarsest
    // level both neighbors intentionally refer to the same single cell.
    int eastIndex(int i) const {return (i+1 == _Nx) ? 0 : i+1;}
    int westIndex(int i) const {return (i == 0) ? _Nx-1 : i-1;}

    double eastBoundaryGhostScale(double vx_face) const
    {
      const double denom = _one_dx - 0.5*vx_face;
      return (_one_dx + 0.5*vx_face) / denom;
    }

    double westBoundaryGhostScale(double vx_face) const
    {
      const double denom = _one_dx + 0.5*vx_face;
      return (_one_dx - 0.5*vx_face) / denom;
    }

    double northBoundaryGhostScale(double vy_face) const
    {
      const double denom = _one_dy - 0.5*vy_face;
      return (_one_dy + 0.5*vy_face) / denom;
    }

    double southBoundaryGhostScale(double vy_face) const
    {
      const double denom = _one_dy + 0.5*vy_face;
      return (_one_dy - 0.5*vy_face) / denom;
    }

    double computeGridMaxNorm(const MatrixXd &v) const
    {
      return v.max_norm();
    }

    void makeProlongation(const MatrixXd &v, MatrixXd &vf) const
    {
      vf.fill(0.0);

      for(int i = 0, i1 = 0; i < _Nx; i++, i1 += 2){
        for(int j = 0, j1 = 0; j < _Ny; j++, j1 += 2){
          double tmp = v[i][j];
          vf[i1][j1] = tmp;
          vf[i1+1][j1] = tmp;
          vf[i1][j1+1] = tmp;
          vf[i1+1][j1+1] = tmp;
        }
      }
    }

    void makeRestriction(const MatrixXd &v, MatrixXd &vc) const
    {
      for(int i = 0, i1 = 0; i < _Nx; i1++, i += 2){
        for(int j = 0, j1 = 0; j < _Ny; j1++, j += 2){
          double sum = v[i][j] + v[i+1][j] + v[i][j+1] + v[i+1][j+1];
          vc[i1][j1] = 0.25*sum;
        }
      }
    }

    void extendToGhostCell(const MatrixXd &v, MatrixXd &v_ex) const
    {
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j++){
          v_ex[i+1][j+1] = v[i][j];
        }
      }
      // Homogeneous total flux remains on the two horizontal walls.
      for(int i = 0; i < _Nx; i++){
        v_ex[i+1][0]       = southBoundaryGhostScale(_vy[i][0]) * v[i][0];
        v_ex[i+1][_Ny+1]   = northBoundaryGhostScale(_vy[i][_Ny]) * v[i][_Ny-1];
      }
      // West/east are one periodic pair, not physical boundaries.
      for(int j = 0; j < _Ny; j++){
        v_ex[0][j+1]       = v[_Nx-1][j];
        v_ex[_Nx+1][j+1]   = v[0][j];
      }
    }
    void removeGhostCell(const MatrixXd &v_ex, MatrixXd &v) const
    {
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j++){
          v[i][j] = v_ex[i+1][j+1];
        }
      }
    }

    void computeMatrixVectorProduct(const MatrixXd &v, MatrixXd &b) const
    {
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j++){

          double uc = v[i][j];
          double ue = v[eastIndex(i)][j];
          double uw = v[westIndex(i)][j];
          double un = j == _Ny-1 ? 0.0       : v[i][j+1];
          double us = j == 0     ? 0.0       : v[i][j-1];

          Array<double, 5> tmp(uc, ue, uw, un, us);
          b[i][j] = _fds_coeff[i][j] * tmp;
        }
      }
    }

    void correctRightHandSide(const MatrixXd &interior_mask,
                              const MatrixXd &crc_func,
                              MatrixXd &b) const
    {
      MatrixXd exterior_mask(1.0+(-1.0)*interior_mask),
               tmp(crc_func);

      computeMatrixVectorProduct(crc_func.cwiseProduct(exterior_mask), tmp);
      b -= tmp.cwiseProduct(interior_mask);

      computeMatrixVectorProduct(crc_func.cwiseProduct(interior_mask), tmp);
      b += tmp.cwiseProduct(exterior_mask);
    }

    void getDiagInv(MatrixXd &one_diag) const
    {
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j++){
          one_diag[i][j] = 1.0/_fds_coeff[i][j][0];
        }
      }
    }
    void getFDSCoeff(Matrix<Array<double, 5>> &coe) const
    {
      // c, e, w, s, n
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j++){
          coe[i][j][0] = 2.0*(_one_dx2+_one_dy2)
                        +0.5*(_vx[i+1][j]-_vx[i][j])*_one_dx
                        +0.5*(_vy[i][j+1]-_vy[i][j])*_one_dy
                        +_kappa;

          coe[i][j][1] = -_one_dx2 + 0.5*_vx[i+1][j]*_one_dx;
          coe[i][j][2] = -_one_dx2 - 0.5*_vx[i][j]*_one_dx;
          coe[i][j][3] = -_one_dy2 + 0.5*_vy[i][j+1]*_one_dy;
          coe[i][j][4] = -_one_dy2 - 0.5*_vy[i][j]*_one_dy;

          // No west/east coefficient elimination: those entries multiply
          // the wrapped periodic neighbors selected by eastIndex/westIndex.
          if (j == _Ny-1) {
            const double gamma_n = northBoundaryGhostScale(_vy[i][j+1]);
            coe[i][j][0] += gamma_n * coe[i][j][3];
            coe[i][j][3] = 0.0;
          }
          if (j == 0) {
            const double gamma_s = southBoundaryGhostScale(_vy[i][j]);
            coe[i][j][0] += gamma_s * coe[i][j][4];
            coe[i][j][4] = 0.0;
          }
        }
      }
    }

    void relaxByJacobi(const MatrixXd &b, MatrixXd &v, double omega = 1.0) const
    {
      MatrixXd w(_Nx, _Ny, 0.0);
      computeMatrixVectorProduct(v, w);
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j++){
          v[i][j] += omega * (b[i][j]-w[i][j])*_one_diag[i][j];
        }
      }
    }

    void relaxByGaussSeidel(const MatrixXd &b, MatrixXd &v, double omega = 1.0) const
    {
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j++){

          double uc = v[i][j];
          double ue = v[eastIndex(i)][j];
          double uw = v[westIndex(i)][j];
          double un = j == _Ny-1 ? v[i][j]   : v[i][j+1];
          double us = j == 0     ? v[i][j]   : v[i][j-1];
          Array<double, 5> tmp(uc, ue, uw, un, us);

          double sum = tmp * _fds_coeff[i][j];
          v[i][j] += omega * (b[i][j] - sum) * _one_diag[i][j];
        }
      }
    }

    void relaxByRBGS(const MatrixXd &b, MatrixXd &v, double omega = 1.0) const
    {
      for(int i = 0; i < _Nx; i++){
        for(int j = 0; j < _Ny; j += 2){
          double uc = v[i][j];
          double ue = v[eastIndex(i)][j];
          double uw = v[westIndex(i)][j];
          double un = j == _Ny-1 ? v[i][j]   : v[i][j+1];
          double us = j == 0     ? v[i][j]   : v[i][j-1];

          Array<double, 5> tmp(uc, ue, uw, un, us);
          double sum = tmp * _fds_coeff[i][j];
          v[i][j] += omega * (b[i][j] - sum) * _one_diag[i][j];
        }
      }

      for(int i = 0; i < _Nx; i++){
        for(int j = 1; j < _Ny; j += 2){
          double uc = v[i][j];
          double ue = v[eastIndex(i)][j];
          double uw = v[westIndex(i)][j];
          double un = j == _Ny-1 ? v[i][j]   : v[i][j+1];
          double us = j == 0     ? v[i][j]   : v[i][j-1];
          Array<double, 5> tmp(uc, ue, uw, un, us);
          double sum = tmp * _fds_coeff[i][j];
          v[i][j] += omega * (b[i][j] - sum) * _one_diag[i][j];
        }
      }
    }

    void relaxByBackGaussSeidel(const MatrixXd &b, MatrixXd &v, double omega = 1.0) const
    {
      for(int i = _Nx-1; i >= 0; i--){
        for(int j = _Ny-1; j >= 0; j--){

          double uc = v[i][j];
          double ue = v[eastIndex(i)][j];
          double uw = v[westIndex(i)][j];
          double un = j == _Ny-1 ? v[i][j]   : v[i][j+1];
          double us = j == 0     ? v[i][j]   : v[i][j-1];

          Array<double, 5> tmp(uc, ue, uw, un, us);
          double sum = tmp*_fds_coeff[i][j];
          v[i][j] += omega * (b[i][j] - sum) * _one_diag[i][j];
        }
      }
    }

    void makeVcycle(const MatrixXd &b, MatrixXd &v) const
    {
      if (0 == coarser) {
        return;
      }

      MatrixXd r(_Nx,_Ny,0.0),
               c(_Nx,_Ny,0.0),
               r_2(_Nx/2,_Ny/2,0.0),
               c_2(_Nx/2,_Ny/2,0.0);

	    const int nu0(2), nu1(2);

	    for(int i = 0; i < nu0; i++){
	    	relaxByGaussSeidel(b, v, 1.0);
	    	//relaxByRBGS(b, v, 1.0);
	    }

	    computeMatrixVectorProduct(v, r);
      r = b - r;

	    makeRestriction(r, r_2);

	    coarser->makeVcycle(r_2, c_2);

	    coarser->makeProlongation(c_2, c);
      v += c;

	    for(int i = 0; i < nu1; i++){
	    	relaxByGaussSeidel(b, v, 1.0);
	    	//relaxByBackGaussSeidel(b, v, 1.0);
	    	//relaxByRBGS(b, v, 1.0);
	    }
    }

    void makeFullMultigrid(const MatrixXd &b, MatrixXd &v) const
    {
      if (0 == coarser) {
        return;
      }

      MatrixXd b_2(_Nx/2, _Ny/2, 0.0),
               v_2(_Nx/2, _Ny/2, 0.0),
               r(_Nx, _Ny, 0.0),
               c(_Nx, _Ny, 0.0);

      makeRestriction(b, b_2);

      coarser->makeFullMultigrid(b_2, v_2);

      coarser->makeProlongation(v_2, v);

      computeMatrixVectorProduct(v, r);
      r = b - r;

      makeVcycle(r, c);
      v += c;
    }

    int solveWithGaussSeidel(const MatrixXd &b, MatrixXd &v,
                             double rtol = 1.0e-8,
                             double atol = 1.0e-15,
                             int maxitrn = 10000) const
    {
      // use update value for convergence chech

      MatrixXd w(_Nx, _Ny, 0.0),
               r(_Nx, _Ny, 0.0),
               c(_Nx, _Ny, 0.0);

      computeMatrixVectorProduct(v, w);
      r = b - w;

      c.fill(0.0);
      relaxByGaussSeidel(r, c);

      double res_norm0 = computeGridMaxNorm(c);
      if (res_norm0 < atol) {
        return 0;
      }

      v += c;

      bool done = false;
      int itrn = 0;
      //std::cout << "itrn = " << itrn << ", res_norm = " << res_norm0 << std::endl;

      for(itrn = 1; itrn < maxitrn; itrn++){

        computeMatrixVectorProduct(c, w);
        r -= w;

        c.fill(0.0);
        relaxByGaussSeidel(r, c);

        v += c;

        double res_norm = computeGridMaxNorm(c);
        //std::cout << "itrn = " << itrn << ", res_norm = " << res_norm << std::endl;

        if (res_norm < res_norm0*rtol || res_norm < atol) {
          done = true;
          break;
        }
      }

      if (done) {
        //std::cout << "#GS = " << itrn << std::endl;
      } else {
        std::cout << "warning: reaching maximum iteration number." << std::endl;
      }
      return itrn;
    }

    int solveWithMultigrid(const MatrixXd &b, MatrixXd &v,
                           double rtol = 1.0e-8,
                           double atol = 1.0e-15,
                           int maxitrn = 100) const
    {
      // use update value for convergence chech

      const int mg_type = 1;

      MatrixXd w(_Nx, _Ny, 0.0),
               r(_Nx, _Ny, 0.0),
               c(_Nx, _Ny, 0.0);

      computeMatrixVectorProduct(v, w);
      r = b - w;

      c.fill(0.0);

      if (mg_type == 1) {
        makeFullMultigrid(r, c);
      } else {
        makeVcycle(r, c);
      }

      double res_norm0 = computeGridMaxNorm(c);
      if (res_norm0 < atol) {
        return 0;
      }

      v += c;

      bool done = false;
      int itrn = 0;
      //std::cout << "itrn = " << itrn << ", res_norm = " << res_norm0 << std::endl;

      for(itrn = 1; itrn < maxitrn; itrn++){

        computeMatrixVectorProduct(c, w);
        r -= w;

        c.fill(0.0);

        if (mg_type == 1) {
          makeFullMultigrid(r, c);
        } else {
          makeVcycle(r, c);
        }

        v += c;

        double res_norm = computeGridMaxNorm(c);
        //std::cout << "itrn = " << itrn << ", res_norm = " << res_norm << std::endl;

        if (res_norm < res_norm0*rtol || res_norm < atol) {
          done = true;
          break;
        }
      }

      if (done) {
        //std::cout << "#Multigrid = " << itrn << std::endl;
      } else {
        std::cout << "warning: reaching maximum iteration number." << std::endl;
      }
      return itrn;
    }

};
