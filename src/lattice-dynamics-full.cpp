/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifdef LATTICE_DYNAMICS
#include <Eigen/Eigenvalues>
#include <complex>
#include "lattice-dynamics-full.h"
#include "box.h"
#include "potential-master.h"
#include "alloc2d.h"
#include <iostream>

LatticeDynamicsFull::LatticeDynamicsFull(PotentialMaster& pm, bool ae) : PotentialCallbackMoleculePhi(pm, false) {
  atomicExp = ae;
  logSum = 0;
  unstable = false;
}

LatticeDynamicsFull::~LatticeDynamicsFull() {
  free2D((void**)waveVectors);
  free(wvCount);
  free3D((void***)matrix);
  free2D((void**)evals);
}

void LatticeDynamicsFull::setNumCells(int x, int y, int z) {
  numCells[0] = x;
  numCells[1] = y;
  numCells[2] = z;

  int nCells = numCells[0]*numCells[1]*numCells[2];
  int numMolecules = box.getTotalNumMolecules();
  int nBasis = numMolecules/nCells;
  if (numMolecules != nCells*nBasis) {
    fprintf(stderr, "warning: number of molecules %d must be a multiple of the number of cells %d\n", numMolecules, nCells);
    abort();
  }
}

void LatticeDynamicsFull::setupForWV(int n, double** wv) {

  int numMolecules = box.getTotalNumMolecules();
  int nCells = numCells[0]*numCells[1]*numCells[2];
  int nBasis = numMolecules/nCells;

  wCount = n;
  waveVectors = (double**)malloc2D(wCount, 3, sizeof(double));

  matrix = (std::complex<double>***)malloc3D(wCount, nmap[nBasis], nmap[nBasis], sizeof(std::complex<double>));
  for (int i=0; i<wCount; i++) {
    for (int j=0; j<nmap[nBasis]; j++) {
      for (int k=0; k<nmap[nBasis]; k++) {
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
}

void LatticeDynamicsFull::compute() {
  PotentialCallbackMoleculePhi::reset();

  int nCells = numCells[0]*numCells[1]*numCells[2];
  int numMolecules = box.getTotalNumMolecules();
  int nBasis = numMolecules/nCells;
  const double* bs = box.getBoxSize();
  double cellSize[3];
  for (int k=0; k<3; k++) cellSize[k] = bs[k]/numCells[k];

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
  double kk[3];
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

  matrix = (std::complex<double>***)malloc3D(wCount, nmap[nBasis], nmap[nBasis], sizeof(std::complex<double>));
  for (int i=0; i<wCount; i++) {
    for (int j=0; j<nmap[nBasis]; j++) {
      for (int k=0; k<nmap[nBasis]; k++) {
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

  //logSum = 0;
  //unstable = false;

  vector<PotentialCallback*> pcs;
  pcs.push_back(this);
  potentialMaster.computeAll(pcs);
}

void LatticeDynamicsFull::allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor) {
  PotentialCallbackMoleculePhi::allComputeFinished(uTot, virialTot, f, virialTensor);
  int nCells = numCells[0]*numCells[1]*numCells[2];
  int numMolecules = box.getTotalNumMolecules();
  int nBasis = numMolecules/nCells;
  const double* bs = box.getBoxSize();
  double cellSize[3];
  for (int k=0; k<3; k++) cellSize[k] = bs[k]/numCells[k];

  int N = box.getTotalNumMolecules();
  std::complex<double> ifac[wCount];
  ifac[0] = 1;
  bool ifacDone[wCount];
  double dr[3];
  evals = (double**)malloc2D(wCount, 3*nBasis, sizeof(double));
  for (int jMolecule=0; jMolecule<N; jMolecule++) {
    double *rj = box.getAtomPosition(jMolecule);
    //for (int g=0; g<3; g++){printf(" %f  ", rj[g]);}
    //printf("\n");
    int jCell = jMolecule/nBasis;
    int jBasis = jMolecule - jCell*nBasis;
    if (jBasis == 0) {
      ifacDone[0] = true;
      for (int i=0; i<wCount; i++) ifacDone[i] = false;
    }
    int xjCell = jCell/(numCells[1]*numCells[2]);
    int yjCell = (jCell - xjCell*numCells[1]*numCells[2]) / numCells[2];
    int zjCell = jCell - xjCell*numCells[1]*numCells[2] - yjCell*numCells[2];
    if (xjCell > numCells[0]/2) xjCell -= numCells[0];
    if (yjCell > numCells[1]/2) yjCell -= numCells[1];
    if (zjCell > numCells[2]/2) zjCell -= numCells[2];
    int xyzjCell[3] = {xjCell,yjCell,zjCell}; 
    for (int iMolecule=0; iMolecule<nBasis; iMolecule++) {
      double *ri = box.getAtomPosition(iMolecule);
      for (int l=0; l<3; l++) dr[l] = rj[l]-ri[l];
      box.nearestImage(dr);
      for (int k=0; k<wCount; k++) {
        if (!ifacDone[k] || atomicExp) {
          double kdotr = 0;
          for (int i=0; i<3; i++) {
            if (atomicExp) {
              kdotr += dr[i]*waveVectors[k][i];
            }else{
              kdotr += cellSize[i]*xyzjCell[i]*waveVectors[k][i];
            }
          }
          std::complex<double> exparg(0, -kdotr);
          ifac[k] = exp(exparg);
          ifacDone[k] = true;
        }
        for (int i=nmap[iMolecule]; i<nmap[iMolecule+1]; i++) {
          for (int j=nmap[jMolecule]; j<nmap[jMolecule+1]; j++) {
            int jj = nmap[jBasis] + (j-nmap[jMolecule]);
            matrix[k][i][jj] += moleculePhiTotal[i][j]*ifac[k];
          }
        }
      }
    }
  }

  for (int k=0; k<wCount; k++) {
    Eigen::MatrixXcd m(nmap[nBasis], nmap[nBasis]);
    for (int i=0; i<nmap[nBasis]; i++) {
      for (int j=0; j<nmap[nBasis]; j++) {
        m(i,j) = matrix[k][i][j];
      }
    }
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigenSolver(m);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd>::EigenvalueType eVals = eigenSolver.eigenvalues();
    double klnsum = 0;
    double k2 = pow(waveVectors[k][0], 2.0) + pow(waveVectors[k][1], 2.0) + pow(waveVectors[k][2], 2.0);
    for (int i=0; i<nmap[nBasis]; i++) {
      //printf("%d %d %25.15e\n", k, i, std::real(eVals(i)));
      double ev = std::real(eVals(i));
      evals[k][i] = ev;
      if (k2>0 || i>2) {
        if (ev <= 0) {
          unstable = true;
          return;
        }
        klnsum += log(ev);
        //printf("%d %d %f (%f) %f\n", k, i, ev, std::imag(eVals(i)), klnsum);
      }
      else {
        //printf("%d %d %f (%f)\n", k, i, ev, std::imag(eVals(i)));
      }
    }
    logSum += wvCount[k]*klnsum;
    //printf("%d %d %f %f\n", k, wvCount[k], klnsum, logSum);
  }
}
#endif
