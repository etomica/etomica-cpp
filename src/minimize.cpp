/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "minimize.h"
#include "potential-callback.h"
#include "potential-master.h"
#include "vector.h"
#include "matrix.h"
#include "rotation-matrix.h"

Minimize::Minimize(PotentialMaster& pm, bool fb) : PotentialCallbackMoleculePhi(pm, fb), stepCount(0) {
  takesVirialTensor = flexBox;
  int numMolecules = box.getTotalNumMolecules();
  int n = nmap[numMolecules];
  if (flexBox) n += 3;
  fMolecule = (double*)malloc(n*sizeof(double));
  lastU = 0;
  lastStep = 0;
  sdStep = false;
  maxDR = 0.3;
  maxDtheta = 0.1;
}

Minimize::~Minimize() {
  free(fMolecule);
}

void Minimize::doStep() {
  reset();
  selfPotentialCallbackVec.clear();
  selfPotentialCallbackVec.push_back(this);
  potentialMaster.computeAll(selfPotentialCallbackVec);
  stepCount++;

  if (sdStep && lastDU>0) lastStep *= 0.2;

  // F, phi
  // F + dF = 0
  // dF = -F
  // dF = - phi dR
  // dR = phi^-1 F
  int N = box.getTotalNumMolecules();
  int nm = N;
  double lastF = 0;
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    double* iF = fMolecule+nmap[iMolecule];
    //printf("F %d %f %f %f  %f %f %f\n", iMolecule, iF[0], iF[1], iF[2], iF[3], iF[4], iF[5]);
    lastF += iF[0]*iF[0] + iF[1]*iF[1] + iF[2]*iF[2] + iF[3]*iF[3] + iF[4]*iF[4] + iF[5]*iF[5];
  }
  lastF = sqrt(lastF)/box.getNumAtoms();
  //printf("RMS F %e\n", lastF);

  Matrix phiMat(nmap[nm], nmap[nm], moleculePhiTotal);
#ifndef MIN_SKIP_UCHECK
  Matrix phiMat0(nmap[nm], nmap[nm]);
  phiMat0.E(phiMat);
  double* fMolecule0 = (double*)malloc(nmap[nm]*sizeof(double));
  for (int i=0; i<nmap[nm]; i++) fMolecule0[i] = fMolecule[i];
#endif
  // We can't solve all 3N translation DOF based on phi and dFdV because
  // the COM mode is free.  Instead enforce that molecule 0 translates to
  // satisfy fixed COM.
  for (int k=0; k<3; k++) {
    for (int j=0; j<N; j++) {
      for (int c=nmap[j]; c<nmap[j+1]; c++) {
        moleculePhiTotal[k][c] = k==(c-nmap[j]) ? 1 : 0;
      }
    }
    fMolecule[k] = 0;
  }

  phiMat.invert();

  phiMat.transform(fMolecule);
  double mr = 0, mt = 0;
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    for (int k=nmap[iMolecule]; k<nmap[iMolecule]+3; k++) {
      if (fabs(fMolecule[k]) > mr) mr = fabs(fMolecule[k]);
    }
    for (int k=nmap[iMolecule]+3; k<nmap[iMolecule+1]; k++) {
      if (fabs(fMolecule[k]) > mt) mt = fabs(fMolecule[k]);
    }
  }
  double scaleback = 1;
  if (mr > maxDR) scaleback = maxDR/mr;
  if (mt > maxDtheta && maxDtheta/mt<scaleback) scaleback = maxDtheta/mt;
  if (scaleback<1) {
    //printf("scaleback %f\n", scaleback);
    for (int i=0; i<nmap[nm]; i++) {
      fMolecule[i] *= scaleback;
    }
  }
  //for  (int i=0; i<nmap[N]; i++) fMolecule[i] = 0;
  //fMolecule[3] = 0.00001;
  //fMolecule[3] = fMolecule[4] = fMolecule[5] = 0;
  //fMolecule[2] = 0.00001; fMolecule[0]=fMolecule[1] = 0;
  //for (int k=0; k<3; k++) fMolecule[k] *= 0.1;
  //printf("dtheta %d %f %f %f\n", 0, fMolecule[0], fMolecule[1], fMolecule[2]);
  sdStep = false;
  double estDU0 = 0;
#ifndef MIN_SKIP_UCHECK
  for (int i=0; i<nmap[nm]; i++) {
    estDU0 -= fMolecule0[i]*fMolecule[i];
    for (int j=0; j<nmap[nm]; j++) {
      estDU0 += 0.5*phiMat0.matrix[i][j]*fMolecule[i]*fMolecule[j];
    }
  }
  //printf("estimated DU: %e\n", estDU/nm);
  double estDU1 = estDU0;
  if (estDU0 > 0) {
    sdStep = true;
    double fmag = 0;
    for (int i=0; i<nmap[nm]; i++) {
      fmag += fMolecule0[i]*fMolecule0[i];
    }
    fmag = sqrt(fmag);
    if (lastStep==0) lastStep = 0.1;
    
    while (estDU1>0) {
      for (int i=0; i<nmap[nm]; i++) {
        fMolecule[i] = fMolecule0[i]/fmag*lastStep;
      }
      estDU1 = 0;
      for (int i=0; i<nmap[nm]; i++) {
        estDU1 -= fMolecule0[i]*fMolecule[i];
        for (int j=0; j<nmap[nm]; j++) {
          estDU1 += 0.5*phiMat0.matrix[i][j]*fMolecule[i]*fMolecule[j];
        }
      }
      //printf("reestimated DU: %e\n", estDU1/nm);
      if (estDU1>0) lastStep *= 0.2;
    }
  }
#endif
  lastDR = 0;
  // fMolecule is now dR
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    double* idR = fMolecule+nmap[iMolecule];
    //printf("m %d %f %f %f\n", iMolecule, idR[0], idR[1], idR[2]);
    RotationMatrix rot;

    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* center = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);

    if (iLastAtom > iFirstAtom) {
      //printf("mr %d %f %f %f\n", iMolecule, idR[3], idR[4], idR[5]);
      double a2 = 0;
      for (int k=0; k<3; k++) a2 += idR[3+k]*idR[3+k];
      double a1 = sqrt(a2);
      double theta = -a1;
      double axis[3];
      for (int k=0; k<3; k++) axis[k] = a1==0 ? (k==0?1:0) : -idR[3+k]/a1;
      //printf("axis %f %f %f  theta %f\n", axis[0], axis[1], axis[2], theta);
      rot.setAxisAngle(axis, theta);
    }
    //printf("%d  % f % f % f  % 12.10f\n", iMolecule, axis[0], axis[1], axis[2], theta);
    //if (iMolecule==0) theta = 0.00000089;
    //axis[0] = 0.02; axis[1] = -0.095; axis[2] = -sqrt(1-axis[1]*axis[1]-axis[0]*axis[0]);
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      double* ri = box.getAtomPosition(iAtom);
  //printf("r %d  %f %f %f\n", iAtom, ri[0], ri[1], ri[2]);
      //printf("%d  %f %f %f  %f %f %f ", iAtom, ri[0], ri[1], ri[2], center[0], center[1], center[2]);
      double idr2[3] = {0,0,0};
      if (iLastAtom>iFirstAtom) {
        //rot.transformAbout(ri, center, box);
        double idr[3];
        for (int k=0; k<3; k++) idr[k] = ri[k] - center[k];
        box.nearestImage(idr);
        for (int k=0; k<3; k++) idr2[k] = idr[k];
        rot.transform(idr);
        //printf("idr %d %f %f %f\n", iAtom, idr[0]-idr2[0], idr[1]-idr2[1], idr[2]-idr2[2]);
        for (int k=0; k<3; k++) {
          idr2[k] = idr[k] - idr2[k];
          ri[k] = center[k] + idr[k];
        }
        //printf("r=> %f %f %f\n", ri[0], ri[1], ri[2]);
      }
      //printf(" %f %f %f ", ri[0], ri[1], ri[2]);
      for (int k=0; k<3; k++) {
        idr2[k] += idR[k];
        lastDR += idr2[k]*idr2[k];
        ri[k] += idR[k];
      }
      box.nearestImage(ri);
    }
  }
  lastDR = sqrt(lastDR)/box.getNumAtoms();
  //printf("step %ld %f %f %f\n", stepCount-1, lastDU, lastDR, lastF);
  potentialMaster.init();
}

void Minimize::doSteps(int steps, double uTol, double rTol) {
  for (int i=0; i<steps; i++) {
    doStep();
  }
}

long Minimize::getStepCount() {
  return stepCount;
}

double Minimize::getLastU() {
  return lastU;
}

double Minimize::getLastDU() {
  return lastDU;
}

double Minimize::getLastDR() {
  return lastDR;
}

void Minimize::allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor) {
  lastDU = uTot - lastU;
  //printf("DU %f\n", lastDU);
  lastU = uTot;
  PotentialCallbackMoleculePhi::allComputeFinished(uTot, virialTot, f, virialTensor);

  // compute molecular force, torque
  int numMolecules = box.getTotalNumMolecules();
  for (int i=0; i<numMolecules; i++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(i, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* ri = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    double atomTorque[3] = {0,0,0};
    for (int k=nmap[i]; k<nmap[i+1]; k++) fMolecule[k] = 0;
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      //double* iF = fMolecule;
      //printf("%d => F %d %f %f %f\n", iAtom, i, iF[0], iF[1], iF[2]);
      for (int k=0; k<3; k++) fMolecule[nmap[i]+k] += f[iAtom][k];
      //printf("%d ==> F %d %f %f %f\n", iAtom, i, iF[0], iF[1], iF[2]);
      if (iLastAtom>iFirstAtom) {
        double* riAtom = box.getAtomPosition(iAtom);
        double drAtom[3] = {riAtom[0]-ri[0], riAtom[1]-ri[1], riAtom[2]-ri[2]};
        box.nearestImage(drAtom);
        Vector::cross(drAtom, f[iAtom], atomTorque);
        //if (i==0) printf("%d torque %f %f %f\n", iAtom, atomTorque[0], atomTorque[1], atomTorque[2]);
        for (int j=0; j<3; j++) fMolecule[nmap[i]+3+j] += atomTorque[j];
      }
    }
  }

  if (flexBox) {
    fMolecule[nmap[numMolecules+1]+0] = virialTensor[0];
    fMolecule[nmap[numMolecules+1]+1] = virialTensor[3];
    fMolecule[nmap[numMolecules+1]+2] = virialTensor[5];
  }
}
