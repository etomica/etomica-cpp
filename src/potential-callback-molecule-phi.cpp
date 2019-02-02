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
PotentialCallbackMoleculePhi::PotentialCallbackMoleculePhi(PotentialMaster& pm, bool fb) : box(pm.getBox()), speciesList(pm.getBox().getSpeciesList()), potentialMaster(pm), flexBox(fb) {
  callFinished = true;
  callPair = true;
  takesForces = true;
  takesPhi = true;
  takesDFDV = flexBox;

  int numAtoms = box.getNumAtoms();
  int N = box.getTotalNumMolecules();
  nmap = (int*)malloc((N+1)*sizeof(int));
  int n = 0;
  for (int i=0; i<N; i++) {
    nmap[i] = n;
    n += 3;
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(i, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    if (iFirstAtom == iLastAtom) continue;
    n += 3;
  }
  nmap[N] = n;
  atomDOF = 3*numAtoms;
  if (flexBox) {
    atomDOF += 3;
    n += 3;
  }
  atomPhiTotal = (double**)malloc2D(atomDOF, atomDOF, sizeof(double));
  moleculePhiTotal = (double**)malloc2D(n, n, sizeof(double));
  com = (double**)malloc2D(N, 3, sizeof(double));
}

PotentialCallbackMoleculePhi::~PotentialCallbackMoleculePhi() {
  free2D((void**)atomPhiTotal);
  free2D((void**)moleculePhiTotal);
  free2D((void**)com);
  free(nmap);
}

void PotentialCallbackMoleculePhi::reset() {
  int N = box.getTotalNumMolecules();
  for (int iMolecule=0; iMolecule<N; iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    double* x = speciesList.get(iSpecies)->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    for (int k=0; k<3; k++) com[iMolecule][k] = x[k];
  }
  int m = 3*box.getNumAtoms();
  if (flexBox) m += 3;
  for (int i=0; i<m; i++) {
    for (int j=0; j<m; j++) {
      atomPhiTotal[i][j] = 0;
    }
  }
}

void PotentialCallbackMoleculePhi::pairComputePhi(int iAtom, int jAtom, double phi[3][3]) {
  for (int k=0; k<3; k++) {
    for (int l=0; l<3; l++) {
      atomPhiTotal[3*iAtom+k][3*jAtom+l] += phi[k][l];
      atomPhiTotal[3*jAtom+l][3*iAtom+k] += phi[k][l];
    }
  }
}

void PotentialCallbackMoleculePhi::pairCompute(int iAtom, int jAtom, double* drij, double u, double du, double d2u) {
  // drij = rj-ri
  double r2 = 0;
  for (int k=0; k<3; k++) r2 += drij[k]*drij[k];
  double dfac = (du - d2u) / (r2*r2);
  for (int k=0; k<3; k++) {
    double foo = drij[k]*dfac;
    for (int l=0; l<3; l++) {
      double der2 = foo*drij[l];
      if (k==l) der2 -= du/r2;
      atomPhiTotal[3*iAtom+k][3*jAtom+l] += der2;
      atomPhiTotal[3*jAtom+l][3*iAtom+k] += der2;
    }
    if (flexBox) {
      // handle this here since pair separation depends directly on V
      int na = box.getNumAtoms();
      for (int l=0; l<3; l++) {
        atomPhiTotal[3*iAtom+k][3*na+l] += dfac * drij[l]*drij[l]*drij[k];
        atomPhiTotal[3*jAtom+k][3*na+l] -= dfac * drij[l]*drij[l]*drij[k];
        atomPhiTotal[3*na+l][3*iAtom+k] += dfac * drij[l]*drij[l]*drij[k];
        atomPhiTotal[3*na+l][3*jAtom+k] -= dfac * drij[l]*drij[l]*drij[k];
        if (k==l) {
          atomPhiTotal[3*iAtom+k][3*na+l] -= du * drij[l] / r2;
          atomPhiTotal[3*jAtom+k][3*na+l] += du * drij[l] / r2;
          atomPhiTotal[3*na+l][3*iAtom+k] -= du * drij[l] / r2;
          atomPhiTotal[3*na+l][3*jAtom+k] += du * drij[l] / r2;
        }
      }
    }
  }
  if (flexBox) {
    int na = box.getNumAtoms();
    for (int k=0; k<3; k++) {
      for (int l=0; l<3; l++) {
        atomPhiTotal[3*na+k][3*na+l] -= dfac * drij[l]*drij[l]*drij[k]*drij[k];
      }
      atomPhiTotal[3*na+k][3*na+k] += 2*du*drij[k]*drij[k]/r2;
    }
  }
}

void PotentialCallbackMoleculePhi::computeDFDV(int iAtom, double* idFdV) {
  int na = box.getNumAtoms();
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  for (int k=0; k<3; k++) {
    for (int l=0; l<3; l++) {
      atomPhiTotal[3*iAtom+k][3*na+l] += idFdV[k]*vol;
      atomPhiTotal[3*na+l][3*iAtom+k] += idFdV[k]*vol;
    }
  }
}

void PotentialCallbackMoleculePhi::allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor) {
  // compute self phi
  int numAtoms = box.getNumAtoms();
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    for (int jAtom=0; jAtom<numAtoms; jAtom++) {
      if (jAtom==iAtom) continue;
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          atomPhiTotal[3*iAtom+k][3*iAtom+l] -= atomPhiTotal[3*iAtom+k][3*jAtom+l];
        }
      }
    }
  }
  int numMolecules = box.getTotalNumMolecules();
  int m = nmap[numMolecules];
  if (flexBox) m += 3;
  for (int i=0; i<m; i++) {
    for (int j=0; j<m; j++) {
      moleculePhiTotal[i][j] = 0;
    }
  }

  if (flexBox) {
    int offset = nmap[numMolecules];
    int na = box.getNumAtoms();
    // add molecular correction to atomic dFdL
    for (int iAtom=0; iAtom<na; iAtom++) {
      for (int jMolecule=0; jMolecule<numMolecules; jMolecule++) {
        int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;
        box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom);
        if (jLastAtom==jFirstAtom) continue;
        Species* species = speciesList.get(jSpecies);
        double* rj = species->getMoleculeCOM(box, jFirstAtom, jLastAtom);
        for (int jAtom=jFirstAtom; jAtom<=jLastAtom; jAtom++) {
          double* rjAtom = box.getAtomPosition(jAtom);
          double drAtom[3] = {rjAtom[0]-rj[0], rjAtom[1]-rj[1], rjAtom[2]-rj[2]};
          box.nearestImage(drAtom);
          for (int k=0; k<3; k++) {
            for (int l=0; l<3; l++) {
              atomPhiTotal[3*iAtom+k][3*na+l] -= atomPhiTotal[3*iAtom+k][3*jAtom+l] * drAtom[l];
              atomPhiTotal[3*na+l][3*iAtom+k] -= atomPhiTotal[3*jAtom+l][3*iAtom+k] * drAtom[l];
            }
          }
        }
      }
    }
    // transform atomic dFdL into molecular dFdL
    for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
      int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
      box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
      Species* species = speciesList.get(iSpecies);
      double* ri = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);
      double atomTorque[3] = {0,0,0};
      for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
        for (int k=0; k<3; k++) {
          for (int l=0; l<3; l++) {
            moleculePhiTotal[nmap[iMolecule]+l][offset+k] += atomPhiTotal[3*iAtom+l][3*na+k];
            moleculePhiTotal[offset+k][nmap[iMolecule]+l] += atomPhiTotal[3*na+k][3*iAtom+l];
          }
        }
        if (iLastAtom>iFirstAtom) {
          double* riAtom = box.getAtomPosition(iAtom);
          double drAtom[3] = {riAtom[0]-ri[0], riAtom[1]-ri[1], riAtom[2]-ri[2]};
          box.nearestImage(drAtom);
          for (int k=0; k<3; k++) {
            Vector::cross(drAtom, atomPhiTotal[3*na+k]+(3*iAtom), atomTorque);
            for (int l=0; l<3; l++) {
              moleculePhiTotal[offset+k][nmap[iMolecule]+l+3] += atomTorque[l];
              moleculePhiTotal[nmap[iMolecule]+l+3][offset+k] += atomTorque[l];
            }
          }
        }
      }
    }
    // and d2U/dLdL
    for (int k=0; k<3; k++) {
      for (int l=0; l<3; l++) {
        moleculePhiTotal[nmap[numMolecules]+k][nmap[numMolecules]+l] = atomPhiTotal[3*na+k][3*na+l];
      }
    }
    // and molecular contribution to d2U/dLdL
    for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
      int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
      box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
      if (iFirstAtom == iLastAtom) continue;
      double* ri = speciesList.get(iSpecies)->getMoleculeCOM(box, iFirstAtom, iLastAtom);
      for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
        double* riAtom = box.getAtomPosition(iAtom);
        double driAtom[3] = {riAtom[0]-ri[0], riAtom[1]-ri[1], riAtom[2]-ri[2]};
        box.nearestImage(driAtom);
        for (int k=0; k<3; k++) {
          moleculePhiTotal[nmap[numMolecules]+k][nmap[numMolecules]+k] -= f[iAtom][k]*driAtom[k];

          for (int jMolecule=0; jMolecule<numMolecules; jMolecule++) {
            int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;
            box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom);
            if (jFirstAtom == jLastAtom) continue;
            double* rj = speciesList.get(jSpecies)->getMoleculeCOM(box, jFirstAtom, jLastAtom);
            for (int jAtom=jFirstAtom; jAtom<=jLastAtom; jAtom++) {
              double* rjAtom = box.getAtomPosition(jAtom);
              double drjAtom[3] = {rjAtom[0]-rj[0], rjAtom[1]-rj[1], rjAtom[2]-rj[2]};
              box.nearestImage(drjAtom);
              for (int l=0; l<3; l++) {
                moleculePhiTotal[nmap[numMolecules]+k][nmap[numMolecules]+l] += atomPhiTotal[3*iAtom+k][3*jAtom+l]*driAtom[k]*drjAtom[l];
              }
            }
          }
        }
      }
    }
  }

  // phi for molecules
  Matrix Ri(3,3), Rj(3,3);
  double** Rim = Ri.matrix;
  Rim[0][0] = Rim[1][1] = Rim[2][2] = 0;
  double** Rjm = Rj.matrix;
  Rjm[0][0] = Rjm[1][1] = Rjm[2][2] = 0;
  Matrix phimat(3,3);
  double** phimm = phimat.matrix;
  Matrix tmpmat(3,3);
  for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
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
          // TT
          for (int k=0; k<3; k++) {
            for (int l=0; l<3; l++) {
              phimm[k][l] = atomPhiTotal[3*iAtom+k][3*jAtom+l];
              moleculePhiTotal[nmap[iMolecule]+k][nmap[jMolecule]+l] += atomPhiTotal[3*iAtom+k][3*jAtom+l];
            }
          }
          // TR
          if (jLastAtom>jFirstAtom) {
            tmpmat.E(phimat);
            tmpmat.TE(Rj);
            for (int k=0; k<3; k++) {
              for (int l=0; l<3; l++) {
                moleculePhiTotal[nmap[iMolecule] + k][nmap[jMolecule] + 3 + l] += tmpmat.matrix[k][l];
              }
            }
          }
          // RT
          if (iLastAtom>iFirstAtom) {
            tmpmat.E(Ri);
            tmpmat.TE(phimat);
            for (int k=0; k<3; k++) {
              for (int l=0; l<3; l++) {
                moleculePhiTotal[nmap[iMolecule] + 3 + k][nmap[jMolecule] + l] += tmpmat.matrix[k][l];
              }
            }
            // RR
            if (jLastAtom>jFirstAtom) {
              tmpmat.TE(Rj);
              for (int k=0; k<3; k++) {
                for (int l=0; l<3; l++) {
                  moleculePhiTotal[nmap[iMolecule] + 3 + k][nmap[jMolecule] + 3 + l] += tmpmat.matrix[k][l];
                }
              }
            }
          }
        }
      }
      // intramolecular RR correction
      if (iLastAtom>iFirstAtom) {
        double xdotf = Vector::dot(dri, f[iAtom]);
        for (int k=0; k<3; k++) {
          moleculePhiTotal[nmap[iMolecule] + 3 + k][nmap[iMolecule] + 3 + k] += xdotf;
          for (int l=0; l<3; l++) {
            moleculePhiTotal[nmap[iMolecule] + 3 + k][nmap[iMolecule] + 3 + l] -= 0.5*(dri[k]*f[iAtom][l] + dri[l]*f[iAtom][k]);
          }
        }
      }
    }
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
}

double** PotentialCallbackMoleculePhi::getMoleculePhi() {
  return moleculePhiTotal;
}
