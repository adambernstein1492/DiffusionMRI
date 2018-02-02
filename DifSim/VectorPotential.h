#include <fftw3.h>
#include <math.h>
#include <vector>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>

typedef std::vector<double>           dvec;

using std::ofstream;
using std::ios;

enum bc_t { ODD = -1, PERIODIC = 0, EVEN = 1 };
typedef fftw_r2r_kind_do_not_use_me fftw_kind;

class CVectorPotential {

  fftw_plan _planR[2][2][2];
  fftw_plan _planK[2][2][2];
  fftw_plan _planRE[2][2][2];
  fftw_plan _planKE[2][2][2];
  fftw_plan _planR2K, _planK2R;
  int       _nx, _ny, _nz;
  int       _nxc, _nyc, _nzc;
  double    _Lx, _Ly, _Lz;
  double    _Nx, _Ny, _Nz;
  bc_t      _BCx0, _BCx1;
  bc_t      _BCy0, _BCy1;
  bc_t      _BCz0, _BCz1;
  bc_t      _bcx, _bcy, _bcz;
  dvec      _kx, _ky, _kz;;
  dvec      _Jr, _Bm;
  dvec      _Gx, _Gy, _Gz;

private:
  int L(int i, int j, int k) { return k+_nz*(j+_ny*i); }
  void printVecXY(dvec &D, int k, const char *op) {
    if (strlen(op) > 0) {
      ofstream ncOFile(op, ios::out);
      for (k=0; k<_nz; k++) {
        for (int i=0; i<_nx; i++) {
          for (int j=0; j<_ny; j++) {
            ncOFile  <<  D[L(i,j,k)];
            if (j==_ny-1) break;
            ncOFile << " ";
          }
          ncOFile << std::endl;
        }
      }
    } else {
      for (int i=0; i<_nx; i++) {
        for (int j=0; j<_ny; j++) {
          std::cout <<  D[L(i,j,k)] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
  void Dx(dvec &D, bc_t b0, bc_t b1) {
    for (int k=0; k<_nz; k++) {
      for (int j=0; j<_ny; j++) {
        if (b0 == PERIODIC) {
          D[L(0,j,k)] = D[L(_nxc,j,k)] = 0;
          for (int i=1; i<_nxc; i++) {
            double tmp = D[L(i,j,k)];
            D[L(i,j,k)] = _kx[i]*D[L(_nx-i,j,k)];
            D[L(_nx-i,j,k)] = _kx[_nx-i]*tmp;
          }
        } else if ((b0 == ODD)&&(b1 == ODD)) {
          for (int i=_nx-1; i>0; i--) {
            D[L(i,j,k)] = _kx[i-1]*D[L(i-1,j,k)];
          }
          D[L(0,j,k)] = 0;
        } else if ((b0 == ODD)&&(b1 == EVEN)) {
          for (int i=0; i<_nx; i++) {
            D[L(i,j,k)] = _kx[i]*D[L(i,j,k)];
          }
        } else if ((b0 == EVEN)&&(b1 == EVEN)) {
          for (int i=1; i<_nx; i++) {
            D[L(i-1,j,k)] = -_kx[i]*D[L(i,j,k)];
          }
          D[L(_nx-1,j,k)] = 0;
        } else if ((b0 == EVEN)&&(b1 == ODD)) {
          for (int i=0; i<_nx; i++) {
            D[L(i,j,k)] = -_kx[i]*D[L(i,j,k)];
          }
        }
        D[L(0,j,k)] = 0;
      }
    }
  }
  void Dy(dvec &D, bc_t b0, bc_t b1) {
    for (int k=0; k<_nz; k++) {
      for (int i=0; i<_nx; i++) {
        if (b0 == PERIODIC) {
          D[L(i,0,k)] = D[L(i,_nyc,k)] = 0;
          for (int j=1; j<_nyc; j++) {
            double tmp = D[L(i,j,k)];
            D[L(i,j,k)] = _ky[j]*D[L(i,_ny-j,k)];
            D[L(i,_ny-j,k)] = _ky[_ny-j]*tmp;
          }
        } else if ((b0 == ODD)&&(b1 == ODD)) {
          for (int j=_ny-1; j>0; j--) {
            D[L(i,j,k)] = _ky[j-1]*D[L(i,j-1,k)];
          }
          D[L(i,0,k)] = 0;
        } else if ((b0 == ODD)&&(b1 == EVEN)) {
          for (int j=0; j<_ny; j++) {
            D[L(i,j,k)] = _ky[j]*D[L(i,j,k)];
          }
        } else if ((b0 == EVEN)&&(b1 == EVEN)) {
          for (int j=1; j<_ny; j++) {
            D[L(i,j-1,k)] = -_ky[j]*D[L(i,j,k)];
          }
          D[L(i,_ny-1,k)] = 0;
        } else if ((b0 == EVEN)&&(b1 == ODD)) {
          for (int j=0; j<_ny; j++) {
            D[L(i,j,k)] = -_ky[j]*D[L(i,j,k)];
          }
        }
        D[L(i,0,k)] = 0;
      }
    }
  }
  void Dz(dvec &D, bc_t b0, bc_t b1) {
    for (int j=0; j<_ny; j++) {
      for (int i=0; i<_nx; i++) {
        if (b0 == PERIODIC) {
          D[L(i,j,0)] = D[L(i,j,_nzc)] = 0;
          for (int k=1; k<_nzc; k++) {
            double tmp = D[L(i,j,k)];
            D[L(i,j,k)] = _kz[k]*D[L(i,j,_nz-k)];
            D[L(i,j,_nz-k)] = _kz[_nz-k]*tmp;
          }
        } else if ((b0 == ODD)&&(b1 == ODD)) {
          for (int k=_nz-1; k>0; k--) {
            D[L(i,j,k)] = _kz[k-1]*D[L(i,j,k-1)];
          }
          D[L(i,j,0)] = 0;
        } else if ((b0 == ODD)&&(b1 == EVEN)) {
          for (int k=0; k<_nz; k++) {
            D[L(i,j,k)] = _kz[k]*D[L(i,j,k)];
          }
        } else if ((b0 == EVEN)&&(b1 == EVEN)) {
          for (int k=1; k<_nz; k++) {
            D[L(i,j,k-1)] = -_kz[k]*D[L(i,j,k)];
          }
          D[L(i,j,_nz-1)] = 0;
        } else if ((b0 == EVEN)&&(b1 == ODD)) {
          for (int k=0; k<_nz; k++) {
            D[L(i,j,k)] = -_kz[k]*D[L(i,j,k)];
          }
        }
        D[L(i,j,0)] = 0;
      }
    }
  }
  void SetType(bc_t &b0, bc_t &b1, dvec &kv, int n, double L, double &N,
               fftw_kind r2k[2], fftw_kind k2r[2]) {
    int nc = n/2;
    if ((b0 == PERIODIC)||(b1 == PERIODIC)) {
      b0 = b1 = PERIODIC;
      r2k[0] = r2k[1] = FFTW_R2HC;
      k2r[0] = k2r[1] = FFTW_HC2R;
      for (int i=1; i<nc; i++) {
        kv[i] = -2.0*M_PI*(double)i/L;
        kv[n-i] = 2.0*M_PI*(double)i/L;
      }
      kv[0] = kv[nc] = 0;
      N = (double)n;
      return;
    }
    if ((b0 == ODD)&&(b1 == ODD)) {
      r2k[0] = FFTW_RODFT10;
      k2r[0] = FFTW_RODFT01;
      r2k[1] = FFTW_REDFT10;
      k2r[1] = FFTW_REDFT01;
      for (int i=0; i<n; i++) {
        kv[i] = M_PI*((double)i+1.0)/L;
      }
      N = 2.0*(double)n;
    } else if ((b0 == ODD)&&(b1 == EVEN)) {
      r2k[0] = FFTW_RODFT11;
      k2r[0] = FFTW_RODFT11;
      r2k[1] = FFTW_REDFT11;
      k2r[1] = FFTW_REDFT11;
      for (int i=0; i<n; i++) {
        kv[i] = M_PI*((double)i+0.5)/L;
      }
      N = 2.0*(double)n;
    } else if ((b0 == EVEN)&&(b1 == EVEN)) {
      r2k[0] = FFTW_REDFT10;
      k2r[0] = FFTW_REDFT01;
      r2k[1] = FFTW_RODFT10;
      k2r[1] = FFTW_RODFT01;
      for (int i=0; i<n; i++) {
        kv[i] = M_PI*((double)i)/L;
      }
      N = 2.0*(double)n;
    } else if ((b0 == EVEN)&&(b1 == ODD)) {
      r2k[0] = FFTW_REDFT11;
      k2r[0] = FFTW_REDFT11;
      r2k[1] = FFTW_RODFT11;
      k2r[1] = FFTW_RODFT11;
      for (int i=0; i<n; i++) {
        kv[i] = M_PI*((double)i+0.5)/L;
      }
      N = 2.0*(double)n;
    }
  }

public:
  CVectorPotential(int nx=10, int ny=10, int nz=10, double Lx=2.0*M_PI,
           double Ly=2.0*M_PI, double Lz=2.0*M_PI,
           bc_t BCx0 = PERIODIC, bc_t BCx1 = PERIODIC,
           bc_t BCy0 = PERIODIC, bc_t BCy1 = PERIODIC,
           bc_t BCz0 = PERIODIC, bc_t BCz1 = PERIODIC) throw( std::bad_alloc) :
    _nx(nx), _ny(ny), _nz(nz), _Lx(Lx), _Ly(Ly), _Lz(Lz),
    _BCx0(BCx0), _BCx1(BCx1), _BCy0(BCy0), _BCy1(BCy1),
    _BCz0(BCz0), _BCz1(BCz1) {

    fftw_kind X2K[2], K2X[2], Y2L[2], L2Y[2], Z2M[2], M2Z[2];
    fftw_kind x2k, y2l, z2m, k2x, l2y, m2z;
    _Jr.resize(_nx*_ny*_nz);
    _Bm.resize(_nx*_ny*_nz);
    _kx.resize(_nx);
    _ky.resize(_ny);
    _kz.resize(_nz);
    _nxc = _nx/2;
    _nyc = _ny/2;
    _nzc = _nz/2;

     SetType(_BCx0, _BCx1, _kx, _nx, _Lx, _Nx, X2K, K2X);
     SetType(_BCy0, _BCy1, _ky, _ny, _Ly, _Ny, Y2L, L2Y);
     SetType(_BCz0, _BCz1, _kz, _nz, _Lz, _Nz, Z2M, M2Z);

    _planR2K=fftw_plan_r2r_3d(_nx,_ny,_nz,&_Jr.front(),&_Bm.front(),
                              X2K[0],Y2L[0],Z2M[0],FFTW_ESTIMATE);
    _planK2R=fftw_plan_r2r_3d(_nx,_ny,_nz,&_Bm.front(),&_Jr.front(),
                              K2X[0],L2Y[0],M2Z[0],FFTW_ESTIMATE);
     for (int i = 0; i < 2; i++)
     for (int j = 0; j < 2; j++)
     for (int k = 0; k < 2; k++) {
      _planR[i][j][k]=fftw_plan_r2r_3d(_nx,_ny,_nz,&_Jr.front(),&_Bm.front(),
                                       X2K[i],Y2L[j],Z2M[k],FFTW_ESTIMATE);
      _planK[i][j][k]=fftw_plan_r2r_3d(_nx,_ny,_nz,&_Bm.front(),&_Jr.front(),
                                       K2X[i],L2Y[j],M2Z[k],FFTW_ESTIMATE);
      if (_BCx0 == PERIODIC) { x2k = FFTW_R2HC; k2x = FFTW_HC2R; _bcx = _BCx0; }
      else if (i==0) { x2k = FFTW_REDFT10; k2x = FFTW_REDFT01; _bcx = EVEN; }
      else { x2k = FFTW_RODFT10; k2x = FFTW_RODFT01; _bcx = EVEN; }
      if (_BCy0 == PERIODIC) { y2l = FFTW_R2HC; l2y = FFTW_HC2R; _bcy = _BCy0; }
      else if (j==0) { y2l = FFTW_REDFT10; l2y = FFTW_REDFT01; _bcy = EVEN; }
      else { y2l = FFTW_RODFT10; l2y = FFTW_RODFT01; _bcy = EVEN; }
      if (_BCz0 == PERIODIC) { z2m = FFTW_R2HC; m2z = FFTW_HC2R; _bcz = _BCz0; }
      else if (k==0) { z2m = FFTW_REDFT10; m2z = FFTW_REDFT01; _bcz = EVEN; }
      else { z2m = FFTW_RODFT10; m2z = FFTW_RODFT01; _bcz = EVEN; }
      _planRE[i][j][k]=fftw_plan_r2r_3d(_nx,_ny,_nz,&_Jr.front(),&_Bm.front(),
                                        x2k,y2l,z2m,FFTW_ESTIMATE);
      _planKE[i][j][k]=fftw_plan_r2r_3d(_nx,_ny,_nz,&_Bm.front(),&_Jr.front(),
                                        k2x,l2y,m2z,FFTW_ESTIMATE);
     }
  }
  double *getJ() { return &_Jr.front(); }
  void printBm( int k, const char *op) { printVecXY(_Bm,k,op); }
  int solveA() throw( std::bad_alloc) {
    _Bm.assign(_Jr.size(),0);
    _Gx.assign(_Jr.size(),0);
    _Gy.assign(_Jr.size(),0);
    _Gz.assign(_Jr.size(),0);

    solveAx();
    solveAy();
    solveAz();
    computeGrad();
    return 0;
  }
  int solveAx() throw( std::bad_alloc) {
    double k2;
    dvec tmp;
    dvec Jk;
    dvec Byk, Bzk, Ax;
    dvec Byr, Bzr;

    tmp.assign(_Jr.begin(),_Jr.end());
    Jk.resize(_nx*_ny*_nz);
    Ax.assign(_Jr.size(),0);
    Byk.assign(_Jr.size(),0);
    Bzk.assign(_Jr.size(),0);
    Byr.assign(_Jr.size(),0);
    Bzr.assign(_Jr.size(),0);

    fftw_execute_r2r(_planR2K,&_Jr.front(),&Jk.front());
    for (int k=0; k<_nz; k++) {
      for (int j=0; j<_ny; j++) {
        for (int i=0; i<_nx; i++) {
          k2 = (_kx[i]*_kx[i]+_ky[j]*_ky[j]+_kz[k]*_kz[k]);
          if (k2 == 0)
            Ax[L(i,j,k)] = 0;
          else
            Ax[L(i,j,k)] = Jk[L(i,j,k)]/k2;
          Byk[L(i,j,k)] = Ax[L(i,j,k)];
          Bzk[L(i,j,k)] =-Ax[L(i,j,k)];
        }
      }
    }
    Dy(Bzk,_BCy0,_BCy1);
    Dz(Byk,_BCz0,_BCz1);
    fftw_execute_r2r(_planK[0][0][1],&Byk.front(),&Byr.front());
    fftw_execute_r2r(_planK[0][1][0],&Bzk.front(),&Bzr.front());
    for (int i=0; i<_nx*_ny*_nz; i++)
      _Bm[i] += sqrt(Byr[i]*Byr[i]+Bzr[i]*Bzr[i])/(_Nx*_Ny*_Nz);
    return 0;
  }
  int solveAy() throw( std::bad_alloc) {
    double k2;
    dvec tmp;
    dvec Jk;
    dvec Bxk, Bzk, Ay;
    dvec Bxr, Bzr;

    tmp.assign(_Jr.begin(),_Jr.end());
    Jk.resize(_nx*_ny*_nz);
    Ay.assign(_Jr.size(),0);
    Bxk.assign(_Jr.size(),0);
    Bzk.assign(_Jr.size(),0);
    Bxr.assign(_Jr.size(),0);
    Bzr.assign(_Jr.size(),0);

    fftw_execute_r2r(_planR2K,&_Jr.front(),&Jk.front());
    for (int k=0; k<_nz; k++) {
      for (int j=0; j<_ny; j++) {
        for (int i=0; i<_nx; i++) {
          k2 = (_kx[i]*_kx[i]+_ky[j]*_ky[j]+_kz[k]*_kz[k]);
          if (k2 == 0)
            Ay[L(i,j,k)] = 0;
          else
            Ay[L(i,j,k)] = Jk[L(i,j,k)]/k2;
          Bxk[L(i,j,k)] =-Ay[L(i,j,k)];
          Bzk[L(i,j,k)] = Ay[L(i,j,k)];
        }
      }
    }
    Dz(Bxk,_BCz0,_BCz1);
    Dx(Bzk,_BCx0,_BCx1);
    fftw_execute_r2r(_planK[0][0][1],&Bxk.front(),&Bxr.front());
    fftw_execute_r2r(_planK[1][0][0],&Bzk.front(),&Bzr.front());
    for (int i=0; i<_nx*_ny*_nz; i++)
      _Bm[i] += sqrt(Bxr[i]*Bxr[i]+Bzr[i]*Bzr[i])/(_Nx*_Ny*_Nz);
    return 0;
  }
  int solveAz() throw( std::bad_alloc) {
    double k2;
    dvec tmp;
    dvec Jk;
    dvec Bxk, Byk, Az;
    dvec Bxr, Byr;

    tmp.assign(_Jr.begin(),_Jr.end());
    Jk.resize(_nx*_ny*_nz);
    Az.assign(_Jr.size(),0);
    Bxk.assign(_Jr.size(),0);
    Byk.assign(_Jr.size(),0);
    Bxr.assign(_Jr.size(),0);
    Byr.assign(_Jr.size(),0);

    fftw_execute_r2r(_planR2K,&_Jr.front(),&Jk.front());
    for (int k=0; k<_nz; k++) {
      for (int j=0; j<_ny; j++) {
        for (int i=0; i<_nx; i++) {
          k2 = (_kx[i]*_kx[i]+_ky[j]*_ky[j]+_kz[k]*_kz[k]);
          if (k2 == 0)
            Az[L(i,j,k)] = 0;
          else
            Az[L(i,j,k)] = Jk[L(i,j,k)]/k2;
          Bxk[L(i,j,k)] = Az[L(i,j,k)];
          Byk[L(i,j,k)] =-Az[L(i,j,k)];
        }
      }
    }
    Dy(Bxk,_BCy0,_BCy1);
    Dx(Byk,_BCx0,_BCx1);
    fftw_execute_r2r(_planK[0][1][0],&Bxk.front(),&Bxr.front());
    fftw_execute_r2r(_planK[1][0][0],&Byk.front(),&Byr.front());
    for (int i=0; i<_nx*_ny*_nz; i++)
      _Bm[i] += sqrt(Bxr[i]*Bxr[i]+Byr[i]*Byr[i])/(_Nx*_Ny*_Nz);
    return 0;
  }
  int computeGrad() throw( std::bad_alloc) {
    dvec Bxk, Byk, Bzk;
    dvec Bm;

    Bxk.assign(_Jr.size(),0);
    Byk.assign(_Jr.size(),0);
    Bzk.assign(_Jr.size(),0);
    Bm.assign(_Bm.begin(),_Bm.end());
    fftw_execute_r2r(_planRE[0][0][0],&Bm.front(),&Bxk.front());
    Byk.assign(Bxk.begin(),Bxk.end());
    Bzk.assign(Bxk.begin(),Bxk.end());
    Dx(Bxk,_bcx,_bcx);
    Dy(Byk,_bcy,_bcy);
    Dz(Bzk,_bcz,_bcz);
    fftw_execute_r2r(_planKE[1][0][0],&Bxk.front(),&_Gx.front());
    fftw_execute_r2r(_planKE[0][1][0],&Byk.front(),&_Gy.front());
    fftw_execute_r2r(_planKE[0][0][1],&Bzk.front(),&_Gz.front());

    for (int i=0; i<_nx*_ny*_nz; i++) {
      _Gx[i] /= (_Nx*_Ny*_Nz);
      _Gy[i] /= (_Nx*_Ny*_Nz);
      _Gz[i] /= (_Nx*_Ny*_Nz);
    }
    return 0;
  }
  void getGradAt(int x0, int y0, int z0,
                 double &Gx, double &Gy, double &Gz) {
    int i0, i1, j0, j1, k0, k1;
    i0 = floor((double)x0/_Lx*_nx);
    i0 = (i0<0)?0:i0;
    i0 = (i0>=_nx)?_nx-1:i0;
    i1 = (i0==_nx-1)?0:i0+1;
    j0 = floor((double)y0/_Ly*_ny);
    j0 = (j0<0)?0:j0;
    j0 = (j0>=_ny)?_ny-1:j0;
    j1 = (j0==_ny-1)?0:j0+1;
    k0 = floor((double)z0/_Lz*_nz);
    k0 = (k0<0)?0:k0;
    k0 = (k0>=_nz)?_nz-1:k0;
    k1 = (k0==_nz-1)?0:k0+1;

    Gz = 0;
    Gx = _Gx[L(i0,j0,k0)]+_Gx[L(i0,j0,k1)]+_Gx[L(i0,j1,k0)]+_Gx[L(i0,j1,k1)]+
         _Gx[L(i1,j0,k0)]+_Gx[L(i1,j0,k1)]+_Gx[L(i1,j1,k0)]+_Gx[L(i1,j1,k1)];
    Gy = _Gy[L(i0,j0,k0)]+_Gy[L(i0,j0,k1)]+_Gy[L(i0,j1,k0)]+_Gy[L(i0,j1,k1)]+
         _Gy[L(i1,j0,k0)]+_Gy[L(i1,j0,k1)]+_Gy[L(i1,j1,k0)]+_Gy[L(i1,j1,k1)];
    Gx /= 8.0;
    Gy /= 8.0;
    Gz /= 8.0;
  }
};
