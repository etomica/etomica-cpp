/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>

#include "potential-master.h"
#include "alloc2d.h"
#include "vector.h"
#include "rotation-matrix.h"
#include "matrix.h"

/**
 * Computes energy with mapped averaging for rigid molecular system.  Also handles atoms.
 */
PotentialCallbackMoleculeHMA::PotentialCallbackMoleculeHMA(Box& b, SpeciesList& sl, PotentialMaster* pm, double T, double Ph) : box(b), speciesList(sl), potentialMaster(pm), temperature(T), Pharm(Ph), returnAnh(false), computingLat(false), computingPshift(false) {
  callFinished = true;
  takesForces = true;
  int N = box.getTotalNumMolecules();
  latticePositions = (double**)malloc2D(N, 3, sizeof(double));
  latticeOrientations = (double**)malloc2D(N, 6, sizeof(double));
  // let's pretend we can't just copy the whole thing at once
  for (int i=0; i<N; i++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(i, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* pos = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    std::copy(pos, pos+3, latticePositions[i]);
    double o[6];
    species->getMoleculeOrientation(box, iFirstAtom, o, o+3);
    std::copy(o, o+6, latticeOrientations[i]);
  }
  data = (double*) malloc(6*sizeof(double));

  dRdV = nullptr;
}

void PotentialCallbackMoleculeHMA::findShiftV() {
  computingPshift = true;
  vector<PotentialCallback*> callbacks;
  callbacks.push_back(this);
  callPair = true;
  takesPhi = true;
  takesDFDV = true;
  int numAtoms = box.getNumAtoms();
  // construct these bits for atoms, transform to molecules later
  phiTotal = (double**)malloc2D(3*numAtoms, 3*numAtoms, sizeof(double));
  dFdV = (double**)malloc2D(numAtoms, 3, sizeof(double));
  com = (double**)malloc2D(numAtoms, 3, sizeof(double));
  for (int i=0; i<numAtoms; i++) {
    for (int k=0; k<3; k++) {
      dFdV[i][k] = 0;
      for (int j=0; j<numAtoms; j++) {
        for (int l=0; l<3; l++) {
          phiTotal[3*i + k][3*j + l] = 0;
        }
      }
    }
  }
  int numMolecules = box.getTotalNumMolecules();
  for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    double* x = speciesList.get(iSpecies)->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    for (int k=0; k<3; k++) com[iMolecule][k] = x[k];
  }
  potentialMaster->computeAll(callbacks);
  computingPshift = false;
  callPair = false;
  takesPhi = false;
  free2D((void**)phiTotal);
  free2D((void**)dFdV);
  free2D((void**)com);
}

PotentialCallbackMoleculeHMA::~PotentialCallbackMoleculeHMA() {
  free2D((void**)latticePositions);
  free2D((void**)latticeOrientations);
  free(data);
  free2D((void**)dRdV);
}

void PotentialCallbackMoleculeHMA::pairComputePhi(int iAtom, int jAtom, double phi[3][3]) {
  for (int k=0; k<3; k++) {
    for (int l=0; l<3; l++) {
      phiTotal[3*iAtom+k][3*jAtom+l] += phi[k][l];
      phiTotal[3*jAtom+l][3*iAtom+k] += phi[k][l];
    }
  }
}

void PotentialCallbackMoleculeHMA::computeDFDV(int iAtom, double* idFdV) {
  for (int k=0; k<3; k++) dFdV[iAtom][k] += idFdV[k];
}

void PotentialCallbackMoleculeHMA::pairCompute(int iAtom, int jAtom, double* drij, double u, double du, double d2u) {
  // drij = rj-ri
  double r2 = 0;
  for (int k=0; k<3; k++) r2 += drij[k]*drij[k];
  double dfac = (du - d2u) / (r2*r2);
  for (int k=0; k<3; k++) {
    double foo = drij[k]*dfac;
    for (int l=0; l<3; l++) {
      double der2 = foo*drij[l];
      if (k==l) der2 -= du/r2;
      phiTotal[3*iAtom+k][3*jAtom+l] += der2;
      phiTotal[3*jAtom+l][3*iAtom+k] += der2;
    }
    // handle this here since pair separation depends directly on V
    dFdV[iAtom][k] += d2u * drij[k] / r2;
    dFdV[jAtom][k] -= d2u * drij[k] / r2;
  }
}

void PotentialCallbackMoleculeHMA::setReturnAnharmonic(bool ra) {
  returnAnh = ra;
  if (ra) {
    computingLat = true;
    vector<PotentialCallback*> callbacks;
    callbacks.push_back(this);
    potentialMaster->computeAll(callbacks);
    computingLat = false;
    uLat = data[0];
    pLat = data[1];
  }
}

int PotentialCallbackMoleculeHMA::getNumData() {return 5;}

void PotentialCallbackMoleculeHMA::computeShift(double** f) {
  // compute self phi
  int numAtoms = box.getNumAtoms();
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    for (int jAtom=0; jAtom<numAtoms; jAtom++) {
      if (jAtom==iAtom) continue;
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          phiTotal[3*iAtom+k][3*iAtom+l] -= phiTotal[3*iAtom+k][3*jAtom+l];
        }
      }
    }
  }

  int numMolecules = box.getTotalNumMolecules();
  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    if (iLastAtom==iFirstAtom) continue;
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      double* ri = box.getAtomPosition(iAtom);
      double dri[3];
      for (int k=0; k<3; k++) dri[k] = ri[k]-com[iMolecule][k];
      box.nearestImage(dri);
      // includes i=j
      for (int jMolecule=0; jMolecule<box.getTotalNumMolecules(); jMolecule++) {
        int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;
        box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom);
        if (jLastAtom==jFirstAtom) continue;
        for (int jAtom=jFirstAtom; jAtom<=jLastAtom; jAtom++) {
          double* rj = box.getAtomPosition(jAtom);
          double drj[3];
          for (int k=0; k<3; k++) drj[k] = rj[k]-com[jMolecule][k];
          box.nearestImage(drj);
          for (int k=0; k<3; k++) {
            for (int l=0; l<3; l++) {
              dFdV[iAtom][k] += phiTotal[3*iAtom+k][3*jAtom+l] * drj[l];
            }
          }
        }
      }
    }
  }
  // transform to molecular coordinates
  
  // dFdV for molecules
  double* Fm = (double*)malloc(6*numMolecules*sizeof(double));
  for (int i=0; i<6*numMolecules; i++) Fm[i]=0;
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  // phi for molecules
  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      for (int k=0; k<3; k++) {
        Fm[6*iMolecule + k] += dFdV[iAtom][k] / (3*vol);
      }
      double* ri = box.getAtomPosition(iAtom);
      double dri[3];
      for (int k=0; k<3; k++) dri[k] = ri[k]-com[iMolecule][k];
      box.nearestImage(dri);
      double torque[3];
      Vector::cross(dri, dFdV[iAtom], torque);
      for (int k=0; k<3; k++) {
        Fm[6*iMolecule + 3 + k] += torque[k] / (3*vol);
      }
    }
    //printf("dFdV %2d  % f % f % f  % f % f % f\n", iMolecule, Fm[6*iMolecule], Fm[6*iMolecule+1], Fm[6*iMolecule+2], Fm[6*iMolecule+3], Fm[6*iMolecule+4], Fm[6*iMolecule+5]);
  }
  Matrix phimall = Matrix(6*numMolecules, 6*numMolecules);
  double** phim = phimall.matrix;
  Matrix Ri(3,3), Rj(3,3);
  double** Rim = Ri.matrix;
  Rim[0][0] = Rim[1][1] = Rim[2][2] = 0;
  double** Rjm = Rj.matrix;
  Rjm[0][0] = Rjm[1][1] = Rjm[2][2] = 0;
  Matrix phimat(3,3);
  double** phimm = phimat.matrix;
  Matrix tmpmat(3,3);
  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      double* ri = box.getAtomPosition(iAtom);
      double dri[3];
      for (int k=0; k<3; k++) dri[k] = ri[k]-com[iMolecule][k];
      box.nearestImage(dri);
      Rim[0][1] = dri[2];
      Rim[1][0] = -dri[2];
      Rim[0][2] = -dri[1];
      Rim[2][0] = dri[1];
      Rim[1][2] = dri[0];
      Rim[2][1] = -dri[0];
      Ri.transpose();
      for (int jMolecule=0; jMolecule<box.getTotalNumMolecules(); jMolecule++) {
        int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;
        box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom);
        for (int jAtom=jFirstAtom; jAtom<=jLastAtom; jAtom++) {
          double* rj = box.getAtomPosition(jAtom);
          double drj[3];
          for (int k=0; k<3; k++) drj[k] = rj[k]-com[jMolecule][k];
          box.nearestImage(drj);
          Rjm[0][1] = drj[2];
          Rjm[1][0] = -drj[2];
          Rjm[0][2] = -drj[1];
          Rjm[2][0] = drj[1];
          Rjm[1][2] = drj[0];
          Rjm[2][1] = -drj[0];
          for (int k=0; k<3; k++) {
            for (int l=0; l<3; l++) {
              phimm[k][l] = phiTotal[3*iAtom+k][3*jAtom+l];
              phim[6*iMolecule+k][6*jMolecule+l] += phiTotal[3*iAtom+k][3*jAtom+l];
            }
          }
          tmpmat.E(phimat);
          tmpmat.TE(Rj);
          for (int k=0; k<3; k++) {
            for (int l=0; l<3; l++) {
              phim[6*iMolecule + k][6*jMolecule + 3 + l] += tmpmat.matrix[k][l];
            }
          }
          tmpmat.E(Ri);
          tmpmat.TE(phimat);
          for (int k=0; k<3; k++) {
            for (int l=0; l<3; l++) {
              phim[6*iMolecule + 3 + k][6*jMolecule + l] += tmpmat.matrix[k][l];
            }
          }
          tmpmat.TE(Rj);
          for (int k=0; k<3; k++) {
            for (int l=0; l<3; l++) {
              phim[6*iMolecule + 3 + k][6*jMolecule + 3 + l] += tmpmat.matrix[k][l];
            }
          }
        }
      }
      double xdotf = Vector::dot(dri, f[iAtom]);
      for (int k=0; k<3; k++) {
        phim[6*iMolecule + 3 + k][6*iMolecule + 3 + k] += xdotf;
        for (int l=0; l<3; l++) {
          phim[6*iMolecule + 3 + k][6*iMolecule + 3 + l] -= 0.5*(dri[k]*f[iAtom][l] + dri[l]*f[iAtom][k]);
        }
      }
    }
  }
  // We can't solve all 3N translation DOF based on phi and dFdV because
  // the COM mode is free.  Instead enforce that molecule 0 translates to
  // satisfy fixed COM.
  for (int k=0; k<3; k++) {
    for (int j=0; j<numMolecules; j++) {
      for (int l=0; l<6; l++) {
        phim[k][6*j+l] = k==l ? 1 : 0;
      }
    }
    Fm[k] = 0;
  }
  /*printf("atomic phi00\n");
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      printf(" % f", phiTotal[i][j]);
    }
    printf("\n");
  }*/
  /*for (int i=0; i<46; i++) {
    printf("phi 45 %d", i);
    for (int k=0; k<3; k++) {
      for (int l=0; l<3; l++) {
        printf(" % 15.8f", phim[6*i + k][6*45+l]);
      }
      printf("\n");
    }
    break;
  }
  printf("phi00\n");
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      printf(" % f", phim[i][j]);
    }
    printf("\n");
  }
  printf("phi01\n");
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      printf(" % f", phim[i][6+j]);
    }
    printf("\n");
  }
  printf("phi11\n");
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      printf(" % f", phim[6+i][6+j]);
    }
    printf("\n");
  }*/
  // -dFdV = phimm dRdV
  // dRdV = - phimm^-1 dFdV
  phimall.invert();
  phimall.transform(Fm);
  // Now determine the shift we need to keep COM fixed
  dRdV = (double**)malloc2D(numMolecules, 6, sizeof(double));
  for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
    for (int k=0; k<6; k++) {
      dRdV[iMolecule][k] = Fm[6*iMolecule+k];
    }
  }
  free(Fm);

  // F_i (new force at old lattice site) = sum phi_ij DRj (amount each j should move)
  // DRj = phi_ii^-1 F_i
  // dRi/dV = phi_ii^-1 dFi/dV
}

double** PotentialCallbackMoleculeHMA::getDRDV() {
  return dRdV;
}

void PotentialCallbackMoleculeHMA::allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor) {
  if (computingPshift) {
    computeShift(f);
    return;
  }
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  if (computingLat) {
    data[0] = uTot;
    data[1] = -virialTot/(3*vol);
    return;
  }
  int N = box.getTotalNumMolecules();
  
  double dr0[3], dr[3];
  int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
  box.getMoleculeInfo(0, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
  Species* species = speciesList.get(iSpecies);
  double* ri = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);
  for (int j=0; j<3; j++) {
    dr0[j] = ri[j] - latticePositions[0][j];
  }
  double fdotdrTot = 0, orientationSum = 0;
  double fdrdVTot = 0, fdrdVTot2 = 0;
  double rotationDOF = 0;
  for (int i=0; i<N; i++) {
    box.getMoleculeInfo(i, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    species = speciesList.get(iSpecies);
    double* ri = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    for (int j=0; j<3; j++) {
      dr[j] = ri[j] - latticePositions[i][j] - dr0[j];
    }
    box.nearestImage(dr);
    double fi[3] = {0,0,0};
    double torque[3] = {0,0,0};
    double atomTorque[3] = {0,0,0};
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      for (int j=0; j<3; j++) fi[j] += f[iAtom][j];
      if (iLastAtom>iFirstAtom) {
        double* riAtom = box.getAtomPosition(iAtom);
        double drAtom[3] = {riAtom[0]-ri[0], riAtom[1]-ri[1], riAtom[2]-ri[2]};
        box.nearestImage(drAtom);
        Vector::cross(drAtom, f[iAtom], atomTorque);
        for (int j=0; j<3; j++) torque[j] += atomTorque[j];
      }
    }
    for (int j=0; j<3; j++) {
      fdotdrTot += fi[j]*dr[j];
      if (dRdV) fdrdVTot += fi[j]*dRdV[i][j];
    }

    double o[9];
    species->getMoleculeOrientation(box, iFirstAtom, o, o+3);
    if (o[0]==0 && o[1]==0 && o[2] == 0) continue;
    double* lo = latticeOrientations[i];
    if (o[3]==0 && o[4]==0 && o[5] == 0) {
      rotationDOF += 2;
      // linear molecule
      double dot = o[0]*lo[0] + o[1]*lo[1] + o[2]*lo[2];
      if (fabs(dot) > 1) {
        dot = dot > 0 ? 1 : -1;
      }
      double axis[3];
      Vector::cross(o, lo, axis);
      double dudt = torque[0]*axis[0] + torque[1]*axis[1] + torque[2]*axis[2];
      orientationSum += (1-dot) / sqrt(1-dot*dot) * dudt;
    }
    else {
      rotationDOF += 3;
      // non-linear molecule
      double lo2[3];
      Vector::cross(lo, lo+3, lo2);
      RotationMatrix rotMatLat;
      rotMatLat.setRows(lo, lo+3, lo2);
      double o2[3];
      Vector::cross(o, o+3, o2);
      RotationMatrix rotMat;
      rotMat.setRows(o, o+3, o2);
      rotMat.transpose();
      rotMat.TE(rotMatLat);
      double axisRot[3];
      axisRot[0] = rotMat.matrix[2][1] - rotMat.matrix[1][2];
      axisRot[1] = rotMat.matrix[0][2] - rotMat.matrix[2][0];
      axisRot[2] = rotMat.matrix[1][0] - rotMat.matrix[0][1];
      double axisLength = sqrt(axisRot[0]*axisRot[0] + axisRot[1]*axisRot[1] + axisRot[2]*axisRot[2]);
      if (axisLength == 0) continue;
      // we'll miss any angle > 90; such angles indicate big trouble anyway
      double theta = asin(0.5 * axisLength);
      for (int j=0; j<3; j++) axisRot[j] /= axisLength;
      double dudt = -(torque[0]*axisRot[0] + torque[1]*axisRot[1] + torque[2]*axisRot[2]);
      if  (dudt == 0) continue;
      double denominator =  1 - cos(theta);
      if (denominator == 0) continue;
      orientationSum -= 1.5 * (theta - 0.5*axisLength) / denominator * dudt;

      if (dRdV) {
        fdrdVTot2 += torque[0]*dRdV[i][3] + torque[1]*dRdV[i][4] + torque[2]*dRdV[i][5];
      }
    }
  }
  double uHarm = 0.5*(3*(N-1) + rotationDOF)*temperature;
  double u0 = returnAnh ? (uHarm + uLat) : 0;
  data[0] = uTot - u0;
  double p0 = returnAnh ? (pLat + Pharm) : 0;
  data[1] = N*temperature/vol - virialTot/(3*vol) - p0;
  data[2] = (returnAnh ? -uLat : uHarm) + uTot + 0.5*fdotdrTot + orientationSum;
  double fV = (Pharm/temperature - N/vol)/(3*(N-1) + rotationDOF);
  //printf("%e %e %e\n", -virialTot/(3*vol), fV * (fdotdrTot + 2*orientationSum), fV);
  data[3] = (returnAnh ? -pLat : Pharm) - virialTot/(3*vol) + fV * (fdotdrTot + 2*orientationSum);
  data[4] = (returnAnh ? -pLat : Pharm) - virialTot/(3*vol) + fV * (fdotdrTot + 2*orientationSum) + fdrdVTot + fdrdVTot2;
}

double* PotentialCallbackMoleculeHMA::getData() {
  return data;
}
