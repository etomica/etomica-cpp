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

EwaldFourier::EwaldFourier(const SpeciesList& sl, Box& b) : EwaldBase(sl,b), kCut(0), alpha(0), eta(0) {
}

EwaldFourier::~EwaldFourier() {
}

void EwaldFourier::setCutoff(double kc) {
  kCut = kc;
}

void EwaldFourier::setChargeAlpha(double a) {
  alpha = a;
}

void EwaldFourier::setR6eta(double e) {
  eta = e;
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

void EwaldFourier::computeAllFourier(const bool doForces, const bool doPhi, const bool doDFDV, double &uTot, double &virialTot, double** force, vector<PotentialCallback*>* pairCallbacks) {
  const int numAtoms = box.getNumAtoms();
  double q2sum = 0, sumBij = 0, sumBii = 0;
  int numAtomTypes = speciesList.getNumAtomTypes();
  for (int iType=0; iType<numAtomTypes; iType++) {
    int inum = numAtomsByType[iType];
    if (inum==0) continue;
    double qi = charges[iType];
    q2sum += inum*qi*qi;
    double Bii = B6[iType][iType];
    if (Bii!=0) {
      sumBii += inum*Bii;
      for (int jType=0; jType<numAtomTypes; jType++) {
        sumBij += inum*numAtomsByType[jType]*B6[iType][jType];
      }
    }
  }
  uTot -= alpha/sqrt(M_PI)*q2sum;
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  if (sumBii > 0) {
    double eta3 = eta*eta*eta;
    uTot -= sqrt(M_PI)*M_PI/(6*vol*eta3) * sumBij;
    virialTot += 3*sqrt(M_PI)*M_PI/(6*vol*eta3) * sumBij;
    uTot += 1.0/(12.0*eta3*eta3) * sumBii;
  }

  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    computeFourierIntramolecular(iMolecule, doForces, doPhi, uTot, virialTot, force);
  }

  const double kCut2 = kCut*kCut;
  int kxMax = (int)(0.5*bs[0]/M_PI*kCut);
  int kMax[3] = {kxMax, (int)(0.5*bs[1]/M_PI*kCut), (int)(0.5*bs[2]/M_PI*kCut)};
  // cube instead of sphere, so conservatively big
  int nk[3] = {kMax[0]+1, 2*kMax[1]+1, 2*kMax[2]+1};
  int nktot = ((2*(nk[0]-1)+1)*nk[1]*nk[2]-1)/2;
  sFac.resize(nktot);
  for (int kB=0; kB<=6; kB++) sFacB[kB].resize(nktot);
  fExp.resize(nktot);
  f6Exp.resize(nktot);
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
      if (charges[iType] == 0 && B6[iType][iType] == 0) continue;
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
  double coeff = 4*M_PI/vol;
  double coeffB = 2*sqrt(M_PI)*M_PI/(24.0*vol);
  double fourierSum = 0, fourierSum6 = 0, virialSum = 0, virialSum6 = 0;
  int ik = 0;
  for (int i=0; i<3; i++) {
    kBasis[i] = 2*M_PI/bs[i];
  }
  sFacAtom.resize(numAtoms);
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
        for (int kB=0; kB<=6; kB++) sFacB[kB][ik] = 0;
        // we could skip this as long as box-length, kCut don't change between calls
        fExp[ik] = 2*coeff*expthing/kxyz2;
        double df6 = 0;
        if (eta>0) {
          double kxyz1 = sqrt(kxyz2);
          double halfketa = 0.5*kxyz1*eta;
          double halfketa2 = halfketa*halfketa;
          f6Exp[ik] = -coeffB*kxyz1*kxyz2*(sqrt(M_PI)*erfc(halfketa) + (0.5/halfketa2 - 1.0)/halfketa*exp(-halfketa2));
          df6 = -(3*f6Exp[ik]/kxyz1 - coeffB*kxyz1*kxyz2*(-exp(-halfketa2)*eta
                  + (-1.5/halfketa2 + 1.0)/halfketa2 * exp(-halfketa2) * 0.5*eta
                  - (0.5/halfketa2 - 1.0)/halfketa * exp(-halfketa2) * 2*0.5*halfketa*eta))*kxyz1;
        }
        for (int iAtom=0; iAtom<numAtoms; iAtom++) {
          int iType = box.getAtomType(iAtom);
          double qi = charges[iType];
          double Bii = B6[iType][iType];
          if (qi==0 && Bii == 0) {
            sFacAtom[iAtom] = 0;
            continue;
          }
          complex<double> icontrib = eik[0][iAtom*nk[0]+ikx]
                          * eik[1][iAtom*nk[1]+kMax[1]+iky]
                          * eik[2][iAtom*nk[2]+kMax[2]+ikz];
          sFacAtom[iAtom] = icontrib;
          sFac[ik] += qi*icontrib;
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacB[kB][ik] += b6[iType][kB]*icontrib;
            }
          }
          if (doPhi) {
            double* ri = box.getAtomPosition(iAtom);
            double phiFac = qi*fExp[ik];
            double phiFac6[7] = {0};
            if (eta>0) {
              for (int kB=0; kB<=6; kB++) {
                phiFac6[kB] = f6Exp[ik]*b6[iType][kB];
              }
            }
            double karray[3] = {kx,ky,kz};
            for (int jAtom=0; jAtom<numAtoms; jAtom++) {
              if (jAtom==iAtom) continue;
              int jType = box.getAtomType(jAtom);
              double qj = charges[jType];
              if (qj==0) continue;
              double* rj = box.getAtomPosition(jAtom);
              double dr[3] = {rj[0]-ri[0],rj[1]-ri[1],rj[2]-ri[2]};
              // cos(kr) ?= f(eik.real,ejk.real)
              double jPhiFac = qj*phiFac;
              if (eta>0) {
                for (int kB=0; kB<=6; kB++) {
                  jPhiFac += b6[jType][6-kB]*phiFac6[kB];
                }
              }
              jPhiFac *= cos(kx*dr[0] + ky*dr[1] + kz*dr[2]);
              for (int k=0; k<3; k++) {
                for (int l=0; l<3; l++) {
                  phi[k][l] = jPhiFac*karray[k]*karray[l];
                }
              }
              if (pairCallbacks) {
                for (vector<PotentialCallback*>::iterator it = pairCallbacks->begin(); it!=pairCallbacks->end(); it++) {
                  (*it)->pairComputePhi(iAtom, jAtom, phi);
                }
              }

            }
          }
        }
        double x = (sFac[ik]*conj(sFac[ik])).real();
        fourierSum += fExp[ik] * x;
        //dFS/dV += dfExp/dV * x;
        //             dfExp/dk dk/dV
        //             dfExp/dk dk/dL dL/dV
        //             dfExp/dk (-2k/3V)
        //             (-2 exp/k^2 - 0.5/alpha^2 exp/k) (-2k/3V)
        //             (2 exp/3V) (2/k + 0.5/alpha^2)
        if (alpha>0) virialSum += 2*expthing * (2/kxyz2 + 0.5/(alpha*alpha))*x;
        // f6Exp = coeffB k^3 (pi^.5 erfc(a) + (0.5/a^3 - 1/a) exp(-a2))
        // dFS6/dV += df6Exp/dV * y
        // df6Exp/dV = df6Exp/dk dk/dV
        //           = df6Exp/dk (-k/3V)
        //           = 3 f6Exp / k + k^3 (pi^.5 (-2 exp(-a^2)/pi^.5) eta/2 + (-1.5/a^4 + 1/a^2) exp(-a2) eta/2 + (0.5/a^4 - 1/a^2) exp(-a2) (-2a) eta/2)
        //           = 3 f6Exp/k + k^3 (-exp(-a^2) eta + (-0.75/a^4 + 0.5/a^2) exp(-a2) eta/a - (0.5/a^4 - 1/a^2) exp(-a2) a eta)
        if (eta>0) {
          for (int kB=0; kB<=6; kB++) {
            double y = (sFacB[kB][ik]*conj(sFacB[6-kB][ik])).real();
            fourierSum6 += f6Exp[ik] * y;
            virialSum6 += df6*y;
          }
        }
        if (doForces) {
          for (int iAtom=0; iAtom<numAtoms; iAtom++) {
            int iType = box.getAtomType(iAtom);
            double coeffki = fExp[ik]*charges[iType]*(sFacAtom[iAtom].imag()*sFac[ik].real()
                             -sFacAtom[iAtom].real()*sFac[ik].imag());
            if (eta>0) {
              for (int kB=0; kB<=7; kB++) {
                coeffki += f6Exp[ik]*b6[iType][kB]*(sFacAtom[iAtom].imag()*sFacB[6-kB][ik].real()
                           -sFacAtom[iAtom].real()*sFacB[6-kB][ik].imag());
              }
            }
            force[iAtom][0] += coeffki * kx;
            force[iAtom][1] += coeffki * ky;
            force[iAtom][2] += coeffki * kz;
            if (doDFDV && pairCallbacks) {
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
  uTot += 0.5*fourierSum;
  // extra factor of 2 here from????
  uTot += fourierSum6;
  // dU/dV = -0.5*coeff*fourierSum/V
  // virialTot = dU/dV *3V
  virialTot += -3*0.5* fourierSum + 0.5*coeff*virialSum  -3*fourierSum6 + virialSum6;
}

double EwaldFourier::uTotalFromAtoms() {
  double q2sum = 0, sumBij = 0, sumBii = 0;
  int numAtomTypes = speciesList.getNumAtomTypes();
  for (int iType=0; iType<numAtomTypes; iType++) {
    int inum = numAtomsByType[iType];
    if (inum==0) continue;
    double qi = charges[iType];
    q2sum += inum*qi*qi;
    double Bii = B6[iType][iType];
    if (Bii!=0) {
      sumBii += inum*Bii;
      for (int jType=0; jType<numAtomTypes; jType++) {
        sumBij += inum*numAtomsByType[jType]*B6[iType][jType];
      }
    }
  }
  double uTot = -alpha/sqrt(M_PI)*q2sum;
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  if (sumBii > 0) {
    double eta3 = eta*eta*eta;
    uTot -= sqrt(M_PI)*M_PI/(6*vol*eta3) * sumBij;
    uTot += 1.0/(12.0*eta3*eta3) * sumBii;
  }

  double virialTot = 0;
  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    computeFourierIntramolecular(iMolecule, false, false, uTot, virialTot, nullptr);
  }
  double fourierSum = 0, fourierSum6 = 0;
  for (int ik=0; ik<(int)sFac.size(); ik++) {
    fourierSum += fExp[ik] * (sFac[ik]*conj(sFac[ik])).real();
    if (eta>0) {
      for (int kB=0; kB<=6; kB++) {
        fourierSum6 += f6Exp[ik] * (sFacB[kB][ik]*conj(sFacB[6-kB][ik])).real();
      }
    }
  }
  uTot += 0.5*fourierSum;
  uTot += fourierSum6;
  return uTot;
}

double EwaldFourier::computeFourierAtom(int iAtom, bool oldEnergy) {
  int iType = box.getAtomType(iAtom);
  double qi = charges[iType];
  if (qi==0 && B6[iType][iType]==0) return 0;

  double u = 0;
  const double kCut2 = kCut*kCut;
  const double* bs = box.getBoxSize();
  int kxMax = (int)(0.5*bs[0]/M_PI*kCut);
  int kMax[3] = {kxMax, (int)(0.5*bs[1]/M_PI*kCut), (int)(0.5*bs[2]/M_PI*kCut)};
  // cube instead of sphere, so conservatively big
  int nk[3] = {kMax[0]+1, 2*kMax[1]+1, 2*kMax[2]+1};
  int nktot = ((2*(nk[0]-1)+1)*nk[1]*nk[2]-1)/2;
  dsFacMolecule.resize(nktot);
  if (eta>0) {
    for (int kB=0; kB<=6; kB++) dsFacBMolecule[kB].resize(nktot);
  }

  double q2Sum = qi*qi;
  u -= alpha/sqrt(M_PI)*q2Sum;
  double vol = bs[0]*bs[1]*bs[2];
  if (B6[iType][iType]!=0) {
    double sumBii = B6[iType][iType];
    double sumBij = 0;
    int numAtomTypes = speciesList.getNumAtomTypes();
    for (int jType=0; jType<numAtomTypes; jType++) {
      sumBij += numAtomsByType[jType]*B6[iType][jType];
    }
    double eta3 = eta*eta*eta;
    u -= sqrt(M_PI)*M_PI/(6*vol*eta3) * sumBij;
    u += 1.0/(12.0*eta3*eta3) * sumBii;
  }
  // we save this too... would need to update after accepted move, revert after rejection
  int numAtoms = box.getNumAtoms();
  for (int a=0; a<3; a++) {
    double fac = 2.0*M_PI/bs[a];
    eik[a].resize(numAtoms*nk[a]); // way bigger than we need, but OK
    if (a==0) q2Sum += qi*qi;
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
  double fourierSum = 0, fourierSum6 = 0;
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
        // when we compute old energy, this will come in as 0
        // when we compute new energy, this will come in as old fourier sum
        if (!oldEnergy) {
          complex<double> sFacMinus = sFac[ik] - dsFacMolecule[ik];
          // energy without this atom
          fourierSum -= fExp[ik] * (sFacMinus*conj(sFacMinus)).real();
          dsFacMolecule[ik] = -dsFacMolecule[ik];
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacMinus = sFacB[kB][ik] - dsFacBMolecule[kB][ik];
              // energy without this atom
              fourierSum6 -= f6Exp[ik] * (sFacMinus*conj(sFacMinus)).real();
              dsFacBMolecule[kB][ik] = -dsFacBMolecule[kB][ik];
            }
          }
        }
        complex<double> icontrib = eik[0][iAtom*nk[0]+ikx]
          * eik[1][iAtom*nk[1]+kMax[1]+iky]
          * eik[2][iAtom*nk[2]+kMax[2]+ikz];
        dsFacMolecule[ik] += qi * icontrib;
        if (eta>0) {
          for (int kB=0; kB<=6; kB++) dsFacBMolecule[kB][ik] += b6[iType][kB]*icontrib;
        }
        
        if (oldEnergy) {
          complex<double> sFacMinus = sFac[ik] - dsFacMolecule[ik];
          // energy with this atom present (as we computed it before) and without
          fourierSum += fExp[ik] * ((sFac[ik]*conj(sFac[ik])).real()
                               -(sFacMinus*conj(sFacMinus)).real());
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacMinus = sFacB[kB][ik] - dsFacBMolecule[kB][ik];
              // energy with this atom present (as we computed it before) and without
              fourierSum6 += f6Exp[ik] * ((sFacB[kB][ik]*conj(sFacB[kB][ik])).real()
                  -(sFacMinus*conj(sFacMinus)).real());
            }
          }
        }
        else {
          complex<double> sFacNew = sFac[ik] + dsFacMolecule[ik];
          fourierSum += fExp[ik] * (sFacNew*conj(sFacNew)).real();
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacNew = sFacB[kB][ik] + dsFacBMolecule[kB][ik];
              fourierSum6 += f6Exp[ik] * (sFacNew*conj(sFacNew)).real();
            }
          }
        }
        ik++;
      }
    }
  }
  u += 0.5*fourierSum;
  u += fourierSum6;
  return u;
}

double EwaldFourier::oneMoleculeFourierEnergy(int iMolecule, bool oldEnergy) {
  double u = 0, v = 0;
  computeFourierIntramolecular(iMolecule, false, false, u, v, nullptr);
  int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);

  const double kCut2 = kCut*kCut;
  const double* bs = box.getBoxSize();
  int kxMax = (int)(0.5*bs[0]/M_PI*kCut);
  int kMax[3] = {kxMax, (int)(0.5*bs[1]/M_PI*kCut), (int)(0.5*bs[2]/M_PI*kCut)};
  // cube instead of sphere, so conservatively big
  int nk[3] = {kMax[0]+1, 2*kMax[1]+1, 2*kMax[2]+1};
  int nktot = ((2*(nk[0]-1)+1)*nk[1]*nk[2]-1)/2;
  dsFacMolecule.resize(nktot);
  if (eta>0) {
    for (int kB=0; kB<=6; kB++) dsFacBMolecule[kB].resize(nktot);
  }

  // we save this too... would need to update after accepted move, revert after rejection
  double q2Sum = 0, sumBij = 0, sumBii = 0;
  int numAtoms = box.getNumAtoms();
  int numAtomTypes = speciesList.getNumAtomTypes();
  for (int a=0; a<3; a++) {
    double fac = 2.0*M_PI/bs[a];
    eik[a].resize(numAtoms*nk[a]); // way bigger than we need, but OK
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      int iType = box.getAtomType(iAtom);
      double qi = charges[iType];
      if (a==0) q2Sum += qi*qi;
      double Bii = B6[iType][iType];
      if (Bii!=0) {
        sumBii += Bii;
        for (int jType=0; jType<numAtomTypes; jType++) {
          sumBij += numAtomsByType[jType]*B6[iType][jType];
        }
      }
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
  if (sumBii == 0 && q2Sum == 0) return 0;
  u -= alpha/sqrt(M_PI)*q2Sum;
  double vol = bs[0]*bs[1]*bs[2];
  if (sumBii > 0) {
    double eta3 = eta*eta*eta;
    u -= sqrt(M_PI)*M_PI/(6*vol*eta3) * sumBij;
    u += 1.0/(12.0*eta3*eta3) * sumBii;
  }
  double fourierSum = 0, fourierSum6 = 0;
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
        // when we compute old energy, this will come in as 0
        // when we compute new energy, this will come in as old fourier sum
        if (!oldEnergy) {
          complex<double> sFacMinus = sFac[ik] - dsFacMolecule[ik];
          // energy without this atom
          fourierSum -= fExp[ik] * (sFacMinus*conj(sFacMinus)).real();
          dsFacMolecule[ik] = -dsFacMolecule[ik];
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacMinus = sFacB[kB][ik] - dsFacBMolecule[kB][ik];
              // energy without this atom
              fourierSum6 -= f6Exp[ik] * (sFacMinus*conj(sFacMinus)).real();
              dsFacBMolecule[kB][ik] = -dsFacBMolecule[kB][ik];
            }
          }
        }
        for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
          int iType = box.getAtomType(iAtom);
          double qi = charges[iType];
          if (qi==0) continue;
          complex<double> icontrib = eik[0][iAtom*nk[0]+ikx]
                          * eik[1][iAtom*nk[1]+kMax[1]+iky]
                          * eik[2][iAtom*nk[2]+kMax[2]+ikz];
          dsFacMolecule[ik] += qi * icontrib;
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) dsFacBMolecule[kB][ik] += b6[iType][kB]*icontrib;
          }
        }
        
        if (oldEnergy) {
          complex<double> sFacMinus = sFac[ik] - dsFacMolecule[ik];
          // energy with this atom present (as we computed it before) and without
          fourierSum += fExp[ik] * ((sFac[ik]*conj(sFac[ik])).real()
                               -(sFacMinus*conj(sFacMinus)).real());
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacMinus = sFacB[kB][ik] - dsFacBMolecule[kB][ik];
              // energy with this atom present (as we computed it before) and without
              fourierSum6 += f6Exp[ik] * ((sFacB[kB][ik]*conj(sFacB[kB][ik])).real()
                  -(sFacMinus*conj(sFacMinus)).real());
            }
          }
        }
        else {
          complex<double> sFacNew = sFac[ik] + dsFacMolecule[ik];
          fourierSum += fExp[ik] * (sFacNew*conj(sFacNew)).real();
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacNew = sFacB[kB][ik] + dsFacBMolecule[kB][ik];
              fourierSum6 += f6Exp[ik] * (sFacNew*conj(sFacNew)).real();
            }
          }
        }
        ik++;
      }
    }
  }
  u += 0.5*fourierSum;
  u += fourierSum6;
  return u;
}

void EwaldFourier::processAtomU(int coeff) {
  for (int i=0; i<(int)sFac.size(); i++) {
    // by the time we get to processAtom(+1), we have
    // the difference.  just use that
    if (coeff==1) {
      sFac[i] += dsFacMolecule[i];
      if (eta>0) {
        for (int kB=0; kB<=6; kB++) sFacB[kB][i] += dsFacBMolecule[kB][i];
      }
    }
    dsFacMolecule[i] = 0;
    if (eta>0) {
      for (int kB=0; kB<=6; kB++) dsFacBMolecule[kB][i] = 0;
    }
  }
}

void EwaldFourier::resetAtomDU() {
  fill(dsFacMolecule.begin(), dsFacMolecule.end(), 0);
  if (eta>0) {
    for (int kB=0; kB<=6; kB++) fill(dsFacBMolecule[kB].begin(), dsFacBMolecule[kB].end(), 0);
  }
}
