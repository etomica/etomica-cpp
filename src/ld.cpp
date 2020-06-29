#ifdef LATTICE_DYNAMICS

#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include "ld.h"
#include "potential.h"
#include "alloc2d.h"

LatticeDynamics::LatticeDynamics(double d, Potential* p, bool lrc, int nb, bool ae) {
  status = 0;

  atomicExp = ae;
  density = d;
  potential = p;
  doLRC = lrc;
  rCut = p->getCutoff();
  lrcFac = 1;
  if (doLRC) {
    lrcFac = 5;
    rCut *= lrcFac;
  }

  if (density < 0.5) {
    fprintf(stderr, "that's a strange density you got there\n");
    status = -1;
  }
  
  double minrc = pow(0.5, 1.0/6.0) / pow(density, 1.0/3.0);
  if (rCut <= minrc) {
    fprintf(stderr, "rc must include first nearest neighbor\n");
    status = -1;
  }

  nBasis = nb;
  basis = (double**)malloc2D(nBasis, 3, sizeof(double));
  cellSize[0] = cellSize[1] = cellSize[2] = 1;
  strain[0] = strain[1] = strain[2] = 1;
}

LatticeDynamics::~LatticeDynamics() {
  if (status < 1) return;
  free2D((void**)waveVectors);
  free3D((void***)matrix);
  free(wvCount);
  free2D((void**)basis);
}

void LatticeDynamics::setUnitCell(double x, double y, double z) {
  cellSize[0] = x;
  cellSize[1] = y;
  cellSize[2] = z;
}

void LatticeDynamics::setNumCells(int x, int y, int z) {
  numCells[0] = x;
  numCells[1] = y;
  numCells[2] = z;
}

void LatticeDynamics::setBasis(int i, double x, double y, double z) {
  basis[i][0] = x;
  basis[i][1] = y;
  basis[i][2] = z;
}

void LatticeDynamics::setStrain(double x, double y, double z) {
  strain[0] = x;
  strain[1] = y;
  strain[2] = z;
}

void LatticeDynamics::setupForWV(int n, double** wv) {

  double vol = cellSize[0] * cellSize[1] * cellSize[2];
  double targetVol = nBasis/density;
  double scale = pow(targetVol/vol, 1.0/3.0);
  for (int k=0; k<3; k++) cellSize[k] *= scale;
  for (int i=0; i<nBasis; i++) {
    for (int k=0; k<3; k++) basis[i][k] *= cellSize[k];
  }

  wCount = n;
  waveVectors = (double**)malloc2D(wCount, 3, sizeof(double));

  matrix = (std::complex<double>***)malloc3D(wCount, 3*nBasis, 3*nBasis, sizeof(std::complex<double>));
  for (int i=0; i<wCount; i++) {
    for (int j=0; j<3*nBasis; j++) {
      for (int k=0; k<3*nBasis; k++) {
        matrix[i][j][k] = 0;
      }
    }
  }

  wCount = n;
  wvCount = (int*)malloc(wCount*sizeof(int));

  for (int i=0; i<n; i++) {
    waveVectors[i][0] = wv[i][0];
    waveVectors[i][1] = wv[i][1];
    waveVectors[i][2] = wv[i][2];
    wvCount[i] = 1;
  }

  selfSum = (double***)malloc3D(nBasis, 3, 3, sizeof(double));
  for (int i=0; i<nBasis; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) selfSum[i][j][k] = 0;
    }
  }

  uLat = 0;
  logSum = 0;
  status = 1;
  for (int i=0; i<3; i++) {
    latticeShells[i] = floor(rCut/cellSize[i])+1;
  }
  rNext = rCut;
  xcellNext = -latticeShells[0];
  ycellNext = -latticeShells[1];
  zcellNext = -latticeShells[2];
  wvNext = wCount-1;
  doneXYZ = 0;
  unstable = false;
}

void LatticeDynamics::setup() {

  double vol = cellSize[0] * cellSize[1] * cellSize[2];
  double targetVol = nBasis/density;
  double scale = pow(targetVol/vol, 1.0/3.0);
  for (int k=0; k<3; k++) cellSize[k] *= scale;
  for (int i=0; i<nBasis; i++) {
    for (int k=0; k<3; k++) basis[i][k] *= cellSize[k];
  }

  int kMin[3], kMax[3];
  for (int i=0; i<3; i++) {
    kMin[i] = (-numCells[i] + 1)/2;
    kMax[i] = numCells[i]/2;
  }
  int*** waveVectorIndices = (int***)malloc3D(2*kMax[0]+1, 2*kMax[1]+1, 2*kMax[2]+1, sizeof(int));
  for (int i=0; i<2*kMax[0]+1; i++) {
    for (int j=0; j<2*kMax[1]+1; j++) {
      for (int k=0; k<2*kMax[2]+1; k++) {
        waveVectorIndices[i][j][k] = 0;
      }
    }
  }
  double wvBasis[3] = {2*M_PI,2*M_PI,2*M_PI};
  for (int k=0; k<3; k++) wvBasis[k] /= cellSize[k]*numCells[k];
  bool flip2[3];
  for (int i=0; i<3; i++) flip2[i] = !(numCells[i]%2);
  Eigen::Vector3i kk;
  wCount = 0;
  for (int k0=kMin[0]; k0<=kMax[0]; k0++) {
    kk[0] = k0;
    for (int k1=kMin[1]; k1<=kMax[1]; k1++) {
      kk[1] = k1;
      for (int k2=kMin[2]; k2<=kMax[2]; k2++) {
        kk[2] = k2;
        int wvIdx[3];
        for (int k=0; k<3; k++) wvIdx[k] = kMax[k];
        int flip = 0;
        for (int i=0; i<3; i++) {
          for (int j=0; j<i; j++) {
            if (kk[j] > 0 && (!flip2[j] || kk[j] < kMax[j])) goto donef;
          }
          if (kk[i] < 0) {
            flip = 1;
            break;
          }
        }

donef:  if (flip) {
          for (int i=0; i<3; i++) {
            wvIdx[i] -= kk[i];
            if (wvIdx[i] == 0 && flip2[i]) {
              wvIdx[i] = 2*kMax[i];
            }
          }
        }
        else {
          for (int i=0; i<3; i++) {
            wvIdx[i] += kk[i];
          }
        }

        if (waveVectorIndices[wvIdx[0]][wvIdx[1]][wvIdx[2]] == 0) {
          wCount++;
        }
        waveVectorIndices[wvIdx[0]][wvIdx[1]][wvIdx[2]]++;
      }
    }
  }

  waveVectors = (double**)malloc2D(wCount, 3, sizeof(double));

  matrix = (std::complex<double>***)malloc3D(wCount, 3*nBasis, 3*nBasis, sizeof(std::complex<double>));
  for (int i=0; i<wCount; i++) {
    for (int j=0; j<3*nBasis; j++) {
      for (int k=0; k<3*nBasis; k++) {
        matrix[i][j][k] = 0;
      }
    }
  }

  wvCount = (int*)malloc(wCount*sizeof(int));
  wCount = 0;

  for (int k0=kMin[0]; k0<=kMax[0]; k0++) {
    for (int k1=kMin[1]; k1<=kMax[1]; k1++) {
      for (int k2=kMin[2]; k2<=kMax[2]; k2++) {
        if (waveVectorIndices[k0+kMax[0]][k1+kMax[1]][k2+kMax[2]] > 0) {
          waveVectors[wCount][0] = k0*wvBasis[0];
          waveVectors[wCount][1] = k1*wvBasis[1];
          waveVectors[wCount][2] = k2*wvBasis[2];
          wvCount[wCount] = waveVectorIndices[k0+kMax[0]][k1+kMax[1]][k2+kMax[2]];
          wCount++;
        }
      }
    }
  }

  free3D((void***)waveVectorIndices);

  selfSum = (double***)malloc3D(nBasis, 3, 3, sizeof(double));
  for (int i=0; i<nBasis; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) selfSum[i][j][k] = 0;
    }
  }

  uLat = 0;
  logSum = 0;
  status = 1;
  for (int i=0; i<3; i++) {
    latticeShells[i] = floor(rCut/cellSize[i])+1;
  }
  rNext = rCut;
  xcellNext = -latticeShells[0];
  ycellNext = -latticeShells[1];
  zcellNext = -latticeShells[2];
  wvNext = wCount-1;
  doneXYZ = 0;
  unstable = false;
}

long long LatticeDynamics::countLS() {
  double d_xyz_cell[3];
  double rCut2 = rCut*rCut;

  long long nxyz = 0;
  int xLS = latticeShells[0];
  for (int xcell=-xLS; xcell<=xLS; xcell++) {
    d_xyz_cell[0] = cellSize[0]*xcell;
    double xMinR2 = xcell==0 ? 0 : (abs(xcell)-1)*cellSize[0];
    xMinR2 *= xMinR2;
    // what is the maximum ycell whose minR would not extend beyond rCut
    // xyMinR2 = xMinR2 + yMinR2 < rCut2
    // yMinR2 < rCut2 - xMinR2
    // yMinR < sqrt(rCut2 - xMinR2)
    // yMinCell = sqrt(rCut2 - xMinR2) / cellSize[1] (round up)
    int yLS = floor(sqrt(rCut2 - xMinR2) / cellSize[1])+1;
    for (int ycell=-yLS; ycell<=yLS; ycell++) {
      d_xyz_cell[1] = cellSize[1]*ycell;
      double yMinR2 = ycell==0 ? 0 : (abs(ycell)-1)*cellSize[1];
      yMinR2 *= yMinR2;
      double xyMinR2 = xMinR2 + yMinR2;
      int zLS = floor(sqrt(rCut2 - xyMinR2) / cellSize[2])+1;
      for (int zcell=-zLS; zcell<=zLS; zcell++) {
        d_xyz_cell[2] = cellSize[2]*zcell;
        int nComp = 0;
        for (int ibasis=0; ibasis<nBasis; ibasis++) {
          for (int jbasis=0; jbasis<nBasis; jbasis++) {
            double r2 = 0;
            double d_xyz[3];
            for (int i=0; i<3; i++) {
              double d_xyz_ij = basis[jbasis][i] - basis[ibasis][i];
              d_xyz[i] = d_xyz_ij + d_xyz_cell[i];
              r2 += d_xyz[i]*d_xyz[i];
            }
            if (r2==0 || r2>rCut2) continue;
            nComp++;
            if (doLRC) {
              nxyz += r2 > rCut2/(lrcFac*lrcFac) ? 2 : (1+wCount);
            }
            else {
              nxyz += 1+wCount;
            }
          }
        }
        //printf("%d %d %d %d\n", xcell, ycell, zcell, nComp);
      }
    }
  }
  return nxyz;
}

long long LatticeDynamics::goLS(int nMax) {

  double d_xyz_cell[3];
  double rCut2 = rCut*rCut;

  if (rNext < 0) {
    fprintf(stderr, "called goLS too many times\n");
    return -1;
  }
  int nxyz = 0;
  std::complex<double> ifac[wCount];
  ifac[0] = 1;
  double der2[3][3];
  for (double rMax = rNext; rMax>0; rMax--) {
    // we want all cells whose minimum distnace (cellMin) follows
    // rMax-1 < cellMin < rMax
    // xCellMin = xcell * cellSize[0]
    // xcell < rMax/cellSize[0]
    // if xcell is very small, we can still reach rMax-1 via y or z
    int xLS = floor(rMax/cellSize[0])+1;
    int xcellStart = -xLS > xcellNext ? -xLS : xcellNext;
    for (int xcell=xcellStart; xcell<=xLS; xcell++) {
      d_xyz_cell[0] = cellSize[0]*xcell;
      double xMinR2 = xcell==0 ? 0 : (abs(xcell)-1)*cellSize[0];
      xMinR2 *= xMinR2;
      int yLS = floor(sqrt(rMax*rMax - xMinR2) / cellSize[1])+1;
      int ycellStart = -yLS > ycellNext ? -yLS : ycellNext;
      for (int ycell=ycellStart; ycell<=yLS; ycell++) {
        d_xyz_cell[1] = cellSize[1]*ycell;
        double yMinR2 = ycell==0 ? 0 : (abs(ycell)-1)*cellSize[1];
        yMinR2 *= yMinR2;
        double xyMinR2 = xMinR2 + yMinR2;
        int zMax = floor(sqrt(rMax*rMax - xyMinR2) / cellSize[2])+1;
        int zMin = 0;
        if (rMax>1) {
          double r2tmp = (rMax-1)*(rMax-1) - xyMinR2;
          if (r2tmp > 0) {
            zMin = floor(sqrt(r2tmp) / cellSize[2])+2;
            if (zMin>zMax) continue;
          }
        }
        for (int zcell=-zMax; zcell<=zMax; zcell++) {
          if (zcell < 0 && zcell > -zMin) zcell = +zMin;
          d_xyz_cell[2] = cellSize[2]*zcell;
          bool ifacDone[wCount];
          ifacDone[0] = true;
          for (int i=1; i<wCount; i++) ifacDone[i] = false;
          int nComp = 0;
          for (int ibasis=0; ibasis<nBasis; ibasis++) {
            for (int jbasis=0; jbasis<nBasis; jbasis++) {
              //if (xcell==0&&ycell==0&&zcell==0&&ibasis==0) printf("jBasis %d\n", jbasis);
              double r2 = 0, r20 = 0;
              double d_xyz[3];
              for (int i=0; i<3; i++) {
                double d_xyz_ij = basis[jbasis][i] - basis[ibasis][i];
                double dx = d_xyz_ij + d_xyz_cell[i];
                r20 += dx*dx;
                d_xyz[i] = strain[i]*dx;
                r2 += d_xyz[i]*d_xyz[i];
              }
              if (r20==0 || r20>rCut2) continue;
              nComp++;
              double u=0, du=0, d2u=0;
              potential->u012(r2, u, du, d2u);
              if (ibasis==0) uLat += u;
              double dfac = (du - d2u) / (r2*r2);
              for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                  der2[i][j] = d_xyz[i]*d_xyz[j] * dfac;
                }
                der2[i][i] -= du/r2;
                for (int j=0; j<3; j++) {
                  selfSum[ibasis][i][j] += der2[i][j];
                }
              }

              nxyz++;
              doneXYZ++;
              for (int k=0; k<wCount; k++) {
                if (doLRC && r2 > rCut2/(lrcFac*lrcFac) && k>0) break;
                nxyz++;
                doneXYZ++;
                if (!ifacDone[k] || atomicExp) {
                  double kdotr = 0;
                  for (int i=0; i<3; i++) {
                    if (atomicExp) {
                      kdotr += d_xyz[i]*waveVectors[k][i];
                    } else {
                      kdotr += d_xyz_cell[i]*waveVectors[k][i];
                    }
                  }
                  std::complex<double> exparg(0, -kdotr);
                  ifac[k] = exp(exparg);
                  ifacDone[k] = true;
                }
                for (int i=0; i<3; i++) {
                  matrix[k][ibasis*3+i][jbasis*3+i] += der2[i][i]*ifac[k];
                  for (int j=i+1; j<3; j++) {
                    matrix[k][ibasis*3+i][jbasis*3+j] += der2[i][j]*ifac[k];
                    matrix[k][ibasis*3+j][jbasis*3+i] += der2[i][j]*ifac[k];
                  }
                }
              }
            }
          }
          //printf("%f %d %d %d %d\n", rMax, xcell, ycell, zcell, nComp);
        }
        if (nxyz > nMax && ycell < yLS) {
          rNext = rMax;
          xcellNext = xcell;
          ycellNext = ycell+1;
          return doneXYZ;
        }
      }
      ycellNext = -latticeShells[1];
      if (nxyz > nMax && xcell < xLS) {
        rNext = rMax;
        xcellNext = xcell+1;
        return doneXYZ;
      }
    }
    xcellNext = -latticeShells[0];
    if (rMax>1) {
      rNext = rMax-1;
      if (rNext<0) rNext = 0;
      if (nxyz > nMax) {
        return doneXYZ;
      }
    }
  }
  return -1;
}

int LatticeDynamics::doSelfSum() {
  for (int ibasis=0; ibasis<nBasis; ibasis++) {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        for (int k=0; k<wCount; k++) {
          matrix[k][ibasis*3+i][ibasis*3+j] -= selfSum[ibasis][i][j];
        }
      }
    }
  }
  free3D((void***)selfSum);
  return -1;
}

double** LatticeDynamics::getEigenvalues(int nwv, int& rwv) {
  if (wvNext<0) {
    fprintf(stderr, "called getEigenvalues too many times");
    rwv = 0;
    return nullptr;
  }
  if (nwv > wvNext+1) nwv = wvNext+1;
  double** rv = (double**)malloc2D(nwv, 3*nBasis, sizeof(double));
  for (int k=wvNext; k>=0; k--) {
    Eigen::MatrixXcd m(nBasis*3, nBasis*3);
    for (int i=0; i<nBasis*3; i++) {
      for (int j=0; j<nBasis*3; j++) {
        m(i,j) = matrix[k][i][j];
      }
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigenSolver(m);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd>::EigenvalueType eVals = eigenSolver.eigenvalues();
    for (int i=0; i<3*nBasis; i++) {
      rv[wvNext-k][i] = std::real(eVals(i));
    }
  }
  wvNext -= nwv;
  if (wvNext<0) wvNext = -1;
  rwv = nwv;
  return rv;
}

int LatticeDynamics::goEVD(int nwv) {
  if (wvNext<0) {
    fprintf(stderr, "called goEVD too many times");
    return -1;
  }
  for (int k=wvNext; k>=0 && k>wvNext-nwv; k--) {
    Eigen::MatrixXcd m(nBasis*3, nBasis*3);
    for (int i=0; i<nBasis*3; i++) {
      for (int j=0; j<nBasis*3; j++) {
        m(i,j) = matrix[k][i][j];
      }
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigenSolver(m);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd>::EigenvalueType eVals = eigenSolver.eigenvalues();
    double klnsum = 0;
    double k2 = pow(waveVectors[k][0], 2.0) + pow(waveVectors[k][1], 2.0) + pow(waveVectors[k][2], 2.0);
    for (int i=0; i<3*nBasis; i++) {
      if (k2>0 || i>2) {
        double ev = std::real(eVals(i));
        if (ev <= 0) {
          unstable = true;
          wvNext = -1;
          return -1;
        }
        klnsum += log(ev);
      }
    }
    logSum += wvCount[k]*klnsum;
  }
  wvNext -= nwv;
  if (wvNext<0) wvNext = -1;
  return wvNext;
}

double LatticeDynamics::getU() {
  double u = 0.5*uLat;
  if (!doLRC) return u;
  double s3 = 1.0/(rCut*rCut*rCut);
  double s9 = s3*s3*s3;
  double A = 4*M_PI;
  double upair = 4.0*A*(s3/3.0);
  double uLRC6 = upair*density/2.0;
  upair = 4.0*A*(s9/9.0);
  double uLRC12 = upair*density/2.0;
  u += uLRC12 - uLRC6;
  return u;
}

#endif
