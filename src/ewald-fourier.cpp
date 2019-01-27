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
  dFdL = (double**)malloc2D(3,3,sizeof(double));
}

EwaldFourier::~EwaldFourier() {
  free2D((void**)dFdL);
}

void EwaldFourier::getOptimalAlpha(double s, double& alpha, double& rc, double& kc) {
  double numAtoms = box.getNumAtoms();
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  // based on crude benchmarks for etomica-cpp
  double tauRatio = 16;
  alpha = pow(tauRatio * M_PI*M_PI*M_PI * numAtoms / (vol*vol), 1.0/6.0);
  rc = s/alpha;
  kc = 2 * alpha*s;
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
    double Bii = B6[iType][iType];
    if (qi==0 && Bii==0) continue;
    double* ri = box.getAtomPosition(iAtom);
    for (int jAtom=iAtom+1; jAtom<=iLastAtom; jAtom++) {
      double qj = charges[box.getAtomType(jAtom)];
      double Bij = B6[iType][jType];
      double qiqj = qi*qj;
      if (qiqj==0 && Bij==0) continue;
      double* rj = box.getAtomPosition(jAtom);
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      box.nearestImage(dr);
      double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      double r = sqrt(r2);
      double ec = 0;
      if (alpha>0) {
        ec = erfc(alpha*r);
        uTot -= qiqj*(1-ec)/r;
      }
      double s2, s4, a2, a4;
      if (eta>0) {
        s2 = 1/r2;
        s4 = s2*s2;
        double eta2 = eta*eta;
        a2 = r2/eta2;
        a4 = a2*a2;
        double ureal = -Bij*eta6r*(1+a2+a4/2)*exp(-a2)/(a4*a2);
        double ufull = -Bij*s4*s2;
        uTot += -ufull + ureal;
      }
      if (!rigidMolecules && doForces) {
        //XXX ensure atom pair is bonded
        double du = 0;
        if (alpha>0) du = -qiqj * (twoosqrtpi * exp(-alpha*alpha*r2) *alpha + (ec-1)/r) / r2;
        if (eta>0) {
          double a6 = a4*a2;
          double dureal = Bij*eta6r*(6 + 6*a2 + 3*a4 + a6)*e/a6;
          double dufull = 6*Bij*s2*s4;
          du += -dufull + dureal;
        }
        for (int m=0; m<3; m++) {
          force[iAtom][m] += dr[m]*du;
          force[jAtom][m] -= dr[m]*du;
        }
        virialTot += du*r2;
      }
    }
  }
}

void EwaldFourier::computeAllFourier(const bool doForces, const bool doPhi, const bool doDFDV, const bool doVirialTensor, double &uTot, double &virialTot, double** force, double* virialTensor, vector<PotentialCallback*>* pairCallbacks) {
  const double alpha6 = eta > 0 ? 1.0/eta : 0;
  const int numAtoms = box.getNumAtoms();
  double q2sum = 0, sumBij = 0, sumBii = 0;
  int numAtomTypes = speciesList.getNumAtomTypes();
  for (int iType=0; iType<numAtomTypes; iType++) {
    int inum = numAtomsByType[iType];
    if (inum==0) continue;
    double qi = charges[iType];
    q2sum += inum*qi*qi;
    double Bii = B6[iType][iType];
    if (eta>0 && Bii!=0) {
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
    double alpha63 = alpha6*alpha6*alpha6;
    uTot -= sqrt(M_PI)*M_PI*alpha63/(6*vol) * sumBij;
    virialTot += 3*sqrt(M_PI)*M_PI*alpha63/(6*vol) * sumBij;
    if (doVirialTensor) {
      virialTensor[0] += sqrt(M_PI)*M_PI*alpha63/(6*vol) * sumBij;
      virialTensor[3] += sqrt(M_PI)*M_PI*alpha63/(6*vol) * sumBij;
      virialTensor[5] += sqrt(M_PI)*M_PI*alpha63/(6*vol) * sumBij;
    }
    uTot += alpha63*alpha63/(12.0) * sumBii;
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
      int idx = iAtom*nk[a];
      if (a>0) idx += kMax[a];
      eik[a][idx] = 1;
      if (nk[a]==1) continue;
      int iType = box.getAtomType(iAtom);
      if (charges[iType] == 0 && B6[iType][iType] == 0) continue;
      double* ri = box.getAtomPosition(iAtom);
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
  double coeffB2 = -2.0*sqrt(M_PI)*M_PI*alpha6*alpha6*alpha6/(3.0*vol);
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
        fExp[ik] = coeff*expthing/kxyz2;
        double hdf6dh = 0;
        if (eta>0) {
          double kxyz1 = sqrt(kxyz2);
          double h = kxyz1/(2*alpha6);
          double h2 = h*h;
          double exph2 = exp(-h2);
          f6Exp[ik] = coeffB2*h*h2*(sqrt(M_PI)*erfc(h) + (0.5/h2 - 1)/h*exph2);
          hdf6dh = 3*f6Exp[ik] - 1.5*coeffB2*exph2;
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
            double iphiFac = 2*qi*fExp[ik];
            double iphiFac6[7] = {0};
            if (eta>0) {
              for (int kB=0; kB<=6; kB++) {
                iphiFac6[kB] = f6Exp[ik]*b6[iType][kB];
              }
            }
            double karray[3] = {kx,ky,kz};
            for (int jAtom=0; jAtom<iAtom; jAtom++) {
              int jType = box.getAtomType(jAtom);
              double qj = charges[jType];
              if (qj*qi==0 && iphiFac6[0]*B6[jType][jType] == 0) continue;
              // cos(kr) ?= f(eik.real,ejk.real)
              double ijPhiFac = qj*iphiFac;
              if (eta>0) {
                for (int kB=0; kB<=6; kB++) {
                  ijPhiFac += b6[jType][6-kB]*iphiFac6[kB];
                }
              }
              ijPhiFac *= 2*(sFacAtom[iAtom]*conj(sFacAtom[jAtom])).real();
              for (int k=0; k<3; k++) {
                for (int l=0; l<3; l++) {
                  phi[k][l] = ijPhiFac*karray[k]*karray[l];
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
        double kdfdk = 0;
        double LdfqdL[3] = {0}, Ldf6dL[3] = {0};
        double dfqdV = 0, df6dV = 0;
        if (alpha>0) {
          kdfdk = -(2 + kxyz2/(2*alpha*alpha))*fExp[ik];
          dfqdV = -fExp[ik]/vol - kdfdk/(3*vol);
          virialSum += 3*vol*dfqdV*x;
          if (doVirialTensor) {
            LdfqdL[0] = -fExp[ik] - kdfdk*kx*kx/kxyz2;
            LdfqdL[1] = -fExp[ik] - kdfdk*ky*ky/kxyz2;
            LdfqdL[2] = -fExp[ik] - kdfdk*kz*kz/kxyz2;
            //double bar = -fExp[ik]*x;
            double foo = -kdfdk*x/kxyz2;
            virialTensor[0] += LdfqdL[0]*x;
            virialTensor[1] += ky*kx * foo;
            virialTensor[2] += kz*kx * foo;
            virialTensor[3] += LdfqdL[1]*x;
            virialTensor[4] += kz*ky * foo;
            virialTensor[5] += LdfqdL[2]*x;
          }
        }
        if (eta>0) {
          df6dV = -f6Exp[ik]/vol - hdf6dh/(3*vol);
          for (int kB=0; kB<=6; kB++) {
            double y = (sFacB[kB][ik]*conj(sFacB[6-kB][ik])).real();
            fourierSum6 += f6Exp[ik] * y;
            virialSum6 += 3*vol*df6dV*y;
            if (doVirialTensor) {
              Ldf6dL[0] = -f6Exp[ik] - hdf6dh*kx*kx/kxyz2;
              Ldf6dL[1] = -f6Exp[ik] - hdf6dh*ky*ky/kxyz2;
              Ldf6dL[2] = -f6Exp[ik] - hdf6dh*kz*kz/kxyz2;
              //double bar = -f6Exp[ik] * y;
              double foo = -hdf6dh*y/kxyz2;
              virialTensor[0] += Ldf6dL[0]*y;
              virialTensor[1] += ky*kx * foo;
              virialTensor[2] += kz*kx * foo;
              virialTensor[3] += Ldf6dL[1]*y;
              virialTensor[4] += kz*ky * foo;
              virialTensor[5] += Ldf6dL[2]*y;
            }
          }
        }
        if (doForces) {
          for (int iAtom=0; iAtom<numAtoms; iAtom++) {
            int iType = box.getAtomType(iAtom);
            double coeffki = alpha==0 ? 0 : 2*fExp[ik]*charges[iType]*(sFacAtom[iAtom]*conj(sFac[ik])).imag();
            double coeffki6 = 0;
            if (eta>0) {
              for (int kB=0; kB<=6; kB++) {
                coeffki6 += 2*f6Exp[ik]*b6[iType][kB]*(sFacAtom[iAtom]*conj(sFacB[6-kB][ik])).imag();
              }
            }
            if (doDFDV && pairCallbacks) {
              double dfdv[3] = {0};
              for (int i=0; i<3; i++) dFdL[i][0] = dFdL[i][1] = dFdL[i][2] = 0;
              double z = 0;
              double kxyz1 = sqrt(kxyz2);
              if (alpha>0) {
                z += 3*vol*(coeffki/fExp[ik])*dfqdV - coeffki;
                // first index is L, second is r
                dFdL[0][0] += LdfqdL[0]*(coeffki/fExp[ik])*kx - coeffki*kx*kx/kxyz1;
                dFdL[0][1] += LdfqdL[0]*(coeffki/fExp[ik])*ky;
                dFdL[0][2] += LdfqdL[0]*(coeffki/fExp[ik])*kz;
                dFdL[1][0] += LdfqdL[1]*(coeffki/fExp[ik])*kx;
                dFdL[1][1] += LdfqdL[1]*(coeffki/fExp[ik])*ky - coeffki*ky*ky/kxyz1;
                dFdL[1][2] += LdfqdL[1]*(coeffki/fExp[ik])*kz;
                dFdL[2][0] += LdfqdL[2]*(coeffki/fExp[ik])*kx;
                dFdL[2][1] += LdfqdL[2]*(coeffki/fExp[ik])*ky;
                dFdL[2][2] += LdfqdL[2]*(coeffki/fExp[ik])*kz - coeffki*kz*kz/kxyz1;
              }
              if (eta>0) {
                z += 3*vol*(coeffki6/f6Exp[ik]) * df6dV - coeffki6;
                dFdL[0][0] += Ldf6dL[0]*(coeffki6/f6Exp[ik])*kx - coeffki6*kx*kx/kxyz1;
                dFdL[0][1] += Ldf6dL[0]*(coeffki6/f6Exp[ik])*ky;
                dFdL[0][2] += Ldf6dL[0]*(coeffki6/f6Exp[ik])*kz;
                dFdL[1][0] += Ldf6dL[1]*(coeffki6/f6Exp[ik])*kx;
                dFdL[1][1] += Ldf6dL[1]*(coeffki6/f6Exp[ik])*ky - coeffki6*ky*ky/kxyz1;
                dFdL[1][2] += Ldf6dL[1]*(coeffki6/f6Exp[ik])*kz;
                dFdL[2][0] += Ldf6dL[2]*(coeffki6/f6Exp[ik])*kx;
                dFdL[2][1] += Ldf6dL[2]*(coeffki6/f6Exp[ik])*ky;
                dFdL[2][2] += Ldf6dL[2]*(coeffki6/f6Exp[ik])*kz - coeffki6*kz*kz/kxyz1;
              }
              dfdv[0] = kx*z;
              dfdv[1] = ky*z;
              dfdv[2] = kz*z;
              for (vector<PotentialCallback*>::iterator it = pairCallbacks->begin(); it!=pairCallbacks->end(); it++) {
                (*it)->computeDFDV(iAtom, dfdv);
                (*it)->computeDFDL(iAtom, dFdL);
              }
            }
            coeffki += coeffki6;
            force[iAtom][0] += coeffki * kx;
            force[iAtom][1] += coeffki * ky;
            force[iAtom][2] += coeffki * kz;
          }
        }
        ik++;
      }
    }
  }
  fill(fExp.begin()+ik, fExp.end(), 0);
  fill(sFac.begin()+ik, sFac.end(), 0);
  if (eta>0) {
    for (int kB=0; kB<=6; kB++) {
      fill(sFacB[kB].begin()+ik, sFacB[kB].end(), 0);
    }
  }
  uTot += fourierSum;
  uTot += fourierSum6;
  virialTot += virialSum + virialSum6;
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
    if (eta>0 && Bii!=0) {
      sumBii += inum*Bii;
      for (int jType=0; jType<numAtomTypes; jType++) {
        sumBij += inum*numAtomsByType[jType]*B6[iType][jType];
      }
    }
  }
  double uTot = -alpha/sqrt(M_PI)*q2sum;
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  if (eta>0 && sumBii > 0) {
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
  uTot += fourierSum;
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
  if (eta > 0 && B6[iType][iType]!=0) {
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
              complex<double> sFacMinus2 = sFacB[6-kB][ik] - dsFacBMolecule[6-kB][ik];
              // energy without this atom
              fourierSum6 -= f6Exp[ik] * (sFacMinus*conj(sFacMinus2)).real();
            }
            for (int kB=0; kB<=6; kB++) {
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
              complex<double> sFacMinus2 = sFacB[6-kB][ik] - dsFacBMolecule[6-kB][ik];
              // energy with this atom present (as we computed it before) and without
              fourierSum6 += f6Exp[ik] * ((sFacB[kB][ik]*conj(sFacB[6-kB][ik])).real()
                  -(sFacMinus*conj(sFacMinus2)).real());
            }
          }
        }
        else {
          complex<double> sFacNew = sFac[ik] + dsFacMolecule[ik];
          fourierSum += fExp[ik] * (sFacNew*conj(sFacNew)).real();
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacNew = sFacB[kB][ik] + dsFacBMolecule[kB][ik];
              complex<double> sFacNew2 = sFacB[6-kB][ik] + dsFacBMolecule[6-kB][ik];
              fourierSum6 += f6Exp[ik] * (sFacNew*conj(sFacNew2)).real();
            }
          }
        }
        ik++;
      }
    }
  }
  fill(dsFacMolecule.begin()+ik, dsFacMolecule.end(), 0);
  if (eta>0) {
    for (int kB=0; kB<=6; kB++) {
      fill(dsFacBMolecule[kB].begin()+ik, dsFacBMolecule[kB].end(), 0);
    }
  }
  u += fourierSum;
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
      int idx = iAtom*nk[a];
      if (a>0) idx += kMax[a];
      eik[a][idx] = 1;
      if (nk[a] == 1) continue;
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
      double* ri = box.getAtomPosition(iAtom);
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
  if (eta>0 && sumBii > 0) {
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
              complex<double> sFacMinus2 = sFacB[6-kB][ik] - dsFacBMolecule[6-kB][ik];
              // energy without this atom
              fourierSum6 -= f6Exp[ik] * (sFacMinus*conj(sFacMinus2)).real();
            }
            for (int kB=0; kB<=6; kB++) {
              dsFacBMolecule[kB][ik] = -dsFacBMolecule[kB][ik];
            }
          }
        }
        for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
          int iType = box.getAtomType(iAtom);
          double qi = charges[iType];
          if (qi==0 && B6[iType][iType] == 0) continue;
          complex<double> icontrib = eik[0][iAtom*nk[0]+ikx]
                          * eik[1][iAtom*nk[1]+kMax[1]+iky]
                          * eik[2][iAtom*nk[2]+kMax[2]+ikz];
          dsFacMolecule[ik] += qi * icontrib;
          if (eta>0 && B6[iType][iType] > 0) {
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
              complex<double> sFacMinus2 = sFacB[6-kB][ik] - dsFacBMolecule[6-kB][ik];
              // energy with this atom present (as we computed it before) and without
              fourierSum6 += f6Exp[ik] * ((sFacB[kB][ik]*conj(sFacB[6-kB][ik])).real()
                  -(sFacMinus*conj(sFacMinus2)).real());
            }
          }
        }
        else {
          complex<double> sFacNew = sFac[ik] + dsFacMolecule[ik];
          fourierSum += fExp[ik] * (sFacNew*conj(sFacNew)).real();
          if (eta>0) {
            for (int kB=0; kB<=6; kB++) {
              sFacNew = sFacB[kB][ik] + dsFacBMolecule[kB][ik];
              complex<double> sFacNew2 = sFacB[6-kB][ik] + dsFacBMolecule[6-kB][ik];
              fourierSum6 += f6Exp[ik] * (sFacNew*conj(sFacNew2)).real();
            }
          }
        }
        ik++;
      }
    }
  }
  fill(dsFacMolecule.begin()+ik, dsFacMolecule.end(), 0);
  if (eta>0) {
    for (int kB=0; kB<=6; kB++) {
      fill(dsFacBMolecule[kB].begin()+ik, dsFacBMolecule[kB].end(), 0);
    }
  }
  u += fourierSum;
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
