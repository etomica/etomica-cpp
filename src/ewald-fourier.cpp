/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <algorithm>
#include <stdio.h>
#include "ewald.h"
#include "box.h"
#include "potential-callback.h"
#include "alloc2d.h"
#include "util.h"

EwaldFourier::EwaldFourier(const SpeciesList& sl, Box& b, vector<PotentialCallback*>* pcb) : EwaldBase(sl,b,pcb) {
  sFacAtom = (complex<double>*)malloc(box.getNumAtoms()*sizeof(complex<double>));
  setEwald(0, 0);
}

EwaldFourier::~EwaldFourier() {
  free(sFacAtom);
}

void EwaldFourier::setEwald(double kc, double a) {
  kCut = kc;
  alpha = a;
  const double* bs = box.getBoxSize();
  for (int i=0; i<3; i++) {
    kBasis[i] = 2*M_PI/bs[i];
  }
}

void EwaldFourier::computeFourierIntramolecular(int iMolecule, const bool doForces, const bool doPhi, double &uTot, double &virialTot, double** force) {
  int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
  if (iLastAtom==iFirstAtom) return;
  double twoosqrtpi = 2.0/sqrt(M_PI);
  for (int iAtom=iFirstAtom; iAtom<iLastAtom; iAtom++) {
    double qi = charges[box.getAtomType(iAtom)];
    if (qi==0) continue;
    double* ri = box.getAtomPosition(iAtom);
    for (int jAtom=iAtom+1; jAtom<=iLastAtom; jAtom++) {
      double qj = charges[box.getAtomType(jAtom)];
      if (qj==0) continue;
      double* rj = box.getAtomPosition(jAtom);
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      box.nearestImage(dr);
      double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      double r = sqrt(r2);
      double qiqj = qi*qj;
      double ec = erfc(alpha*r);
      uTot -= qiqj*(1-ec)/r;
      if (!rigidMolecules && doForces) {
        double du = -qiqj * (twoosqrtpi * exp(-alpha*alpha*r2) *alpha + (ec-1)/r) / r2;
        for (int m=0; m<3; m++) {
          force[iAtom][m] += dr[m]*du;
          force[jAtom][m] -= dr[m]*du;
        }
        virialTot += du*r2;
      }
    }
  }
}

void EwaldFourier::computeAllFourier(const bool doForces, const bool doPhi, const bool doDFDV, double &uTot, double &virialTot, double** force) {
  const int numAtoms = box.getNumAtoms();
  double q2sum = 0;
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double qi = charges[iType];
    if (qi==0) continue;
    q2sum += qi*qi;
  }
  uTot -= alpha/sqrt(M_PI)*q2sum;

  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    computeFourierIntramolecular(iMolecule, doForces, doPhi, uTot, virialTot, force);
  }

  const double kCut2 = kCut*kCut;
  const double* bs = box.getBoxSize();
  int kxMax = (int)(0.5*bs[0]/M_PI*kCut);
  int kMax[3] = {kxMax, (int)(0.5*bs[1]/M_PI*kCut), (int)(0.5*bs[2]/M_PI*kCut)};
  // cube instead of sphere, so conservatively big
  int nk[3] = {kMax[0]+1, 2*kMax[1]+1, 2*kMax[2]+1};
  int nktot = ((2*(nk[0]-1)+1)*nk[1]*nk[2]-1)/2;
  sFac.resize(nktot);
  fExp.resize(nktot);
  // We want exp(i dot(k,r)) for every k and every r
  // then sum over atoms, s(k) = sum[exp(i dot(k,r))] and U(k) = s(k) * s*(k)
  //
  // To get there, we first invoke that
  // exp(i dot(k,r)) = exp(i kx rx) exp(i ky ry) exp(i kz rz)
  // the values of ka (a=x,y,z) will be such that ka = 2 l pi / bs[a]
  // so ea(l=2) = ea(l=1) ea(l=1)
  //    ea(l=3) = ea(l=1) ea(l=2)
  // and so on, where ea(l) = exp(i (2 l pi / bs[a]) ra)
  // See Allen & Tildesley for more details
  // https://dx.doi.org/10.1093/oso/9780198803195.001.0001
  // https://github.com/Allen-Tildesley/examples/blob/master/ewald_module.f90
  for (int a=0; a<3; a++) {
    double fac = 2.0*M_PI/bs[a];
    eik[a].resize(numAtoms*nk[a]);
    for (int iAtom=0; iAtom<numAtoms; iAtom++) {
      int iType = box.getAtomType(iAtom);
      if (charges[iType] == 0) continue;
      int idx = iAtom*nk[a];
      if (a>0) idx += kMax[a];
      double* ri = box.getAtomPosition(iAtom);
      eik[a][idx] = 1;
      eik[a][idx+1] = std::complex<double>(cos(fac*ri[a]), sin(fac*ri[a]));
      for (int i=2; i<=kMax[a]; i++) {
        eik[a][idx+i] = eik[a][idx+1] * eik[a][idx+i-1];
      }
      if (a==0) continue;
      for (int i=1; i<=kMax[a]; i++) {
        eik[a][idx-i] = conj(eik[a][idx+i]);
      }
    }
  }
  double vol = bs[0]*bs[1]*bs[2];
  double coeff = 4*M_PI/vol;
  double fourierSum = 0, virialSum = 0;;
  int ik = 0;
  for (int i=0; i<3; i++) {
    kBasis[i] = 2*M_PI/bs[i];
  }
  for (int ikx=0; ikx<=kxMax; ikx++) {
    double kx = ikx*kBasis[0];
    double kx2 = kx*kx;
    double kyCut2 = kCut2 - kx2;
    bool xpositive = ikx>0;
    int kyMax = (int)(0.5*bs[1]*sqrt(kyCut2)/M_PI);
    for (int iky=-kyMax; iky<=kyMax; iky++) {
      if (!xpositive && iky<0) continue;
      bool ypositive = iky>0;
      double ky = iky*kBasis[1];
      double kxy2 = kx2 + ky*ky;
      int kzMax = (int)(0.5*bs[2]*sqrt(kCut2 - kxy2)/M_PI);
      for (int ikz=-kzMax; ikz<=kzMax; ikz++) {
        if (!xpositive && !ypositive && ikz<=0) continue;
        double kz = ikz*kBasis[2];
        double kxyz2 = kxy2 + kz*kz;
        double expthing = exp(-0.25*kxyz2/(alpha*alpha));
        sFac[ik] = 0;
        for (int iAtom=0; iAtom<numAtoms; iAtom++) {
          int iType = box.getAtomType(iAtom);
          double qi = charges[iType];
          if (qi==0) {
            sFacAtom[iAtom] = 0;
            continue;
          }
          sFacAtom[iAtom] = qi * eik[0][iAtom*nk[0]+ikx]
                               * eik[1][iAtom*nk[1]+kMax[1]+iky]
                               * eik[2][iAtom*nk[2]+kMax[2]+ikz];
          sFac[ik] += sFacAtom[iAtom];
          if (doPhi) {
            double* ri = box.getAtomPosition(iAtom);
            double phiFac = 4*M_PI*qi*expthing/(kxyz2*vol);
            double karray[3] = {kx,ky,kz};
            for (int jAtom=0; jAtom<numAtoms; jAtom++) {
              if (jAtom==iAtom) continue;
              int jType = box.getAtomType(jAtom);
              double qj = charges[jType];
              if (qj==0) continue;
              double* rj = box.getAtomPosition(jAtom);
              double dr[3] = {rj[0]-ri[0],rj[1]-ri[1],rj[2]-ri[2]};
              double jPhiFac = qj*phiFac*cos(kx*dr[0] + ky*dr[1] + kz*dr[2]);
              for (int k=0; k<3; k++) {
                for (int l=0; l<3; l++) {
                  phi[k][l] = jPhiFac*karray[k]*karray[l];
                }
              }
              for (vector<PotentialCallback*>::iterator it = pairCallbacks->begin(); it!=pairCallbacks->end(); it++) {
                (*it)->pairComputePhi(iAtom, jAtom, phi);
              }

            }
          }
        }
        // we could skip this as long as box-length, kCut don't change between calls
        fExp[ik] = 2*expthing/kxyz2;
        double x = (sFac[ik]*conj(sFac[ik])).real();
        fourierSum += fExp[ik] * x;
        //dFS/dV += dfExp/dV * x;
        //             dfExp/dk dk/dV
        //             dfExp/dk dk/dL dL/dV
        //             dfExp/dk (-2k/3V)
        //             (-2 exp/k^2 - 0.5/alpha^2 exp/k) (-2k/3V)
        //             (2 exp/3V) (2/k + 0.5/alpha^2)
        virialSum += 2*expthing * (2/kxyz2 + 0.5/(alpha*alpha))*x;
        if (doForces) {
          double coeffk = coeff * fExp[ik];
          for (int iAtom=0; iAtom<numAtoms; iAtom++) {
            double coeffki = coeffk * (sFacAtom[iAtom].imag()*sFac[ik].real()
                                      -sFacAtom[iAtom].real()*sFac[ik].imag());
            force[iAtom][0] += coeffki * kx;
            force[iAtom][1] += coeffki * ky;
            force[iAtom][2] += coeffki * kz;
            if (doDFDV) {
              double dfdv[3] = {coeffki*kx/(3*vol), coeffki*ky/(3*vol), coeffki*kz/(3*vol)};
              for (vector<PotentialCallback*>::iterator it = pairCallbacks->begin(); it!=pairCallbacks->end(); it++) {
                (*it)->computeDFDV(iAtom, dfdv);
              }
            }
          }
        }
        ik++;
      }
    }
  }
  fill(fExp.begin()+ik, fExp.end(), 0);
  fill(sFac.begin()+ik, sFac.end(), 0);
  uTot += 0.5*coeff * fourierSum;
  // dU/dV = -0.5*coeff*fourierSum/V
  // virialTot = dU/dV *3V
  virialTot += -3*0.5*coeff * fourierSum + 0.5*coeff*virialSum;
}

double EwaldFourier::uTotalFromAtoms() {
  const int numAtoms = box.getNumAtoms();
  double q2sum = 0;
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double qi = charges[iType];
    if (qi==0) continue;
    q2sum += qi*qi;
  }
  double uTot = alpha/sqrt(M_PI)*q2sum;

  double virialTot = 0;
  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    computeFourierIntramolecular(iMolecule, false, false, uTot, virialTot, nullptr);
  }
  double fourierSum = 0;
  for (int ik=0; ik<(int)sFac.size(); ik++) {
    fourierSum += fExp[ik] * (sFac[ik]*conj(sFac[ik])).real();
  }
  const double* bs = box.getBoxSize();
  double coeff = 4*M_PI/(bs[0]*bs[1]*bs[2]);
  uTot += 0.5*coeff * fourierSum;
  return uTot;
}
