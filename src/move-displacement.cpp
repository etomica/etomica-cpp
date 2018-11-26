/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move.h"

MCMoveDisplacement::MCMoveDisplacement(Box& b, PotentialMaster& p, Random& r, double ss) : MCMove(b,p,r,ss) {
}

MCMoveDisplacement::~MCMoveDisplacement() {}

bool MCMoveDisplacement::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int na = box.getNumAtoms();
  if (na==0) {
    iAtom = -1;
    return false;
  }
  iAtom = random.nextInt(na);
  uOld = potentialMaster.oldEnergy(iAtom);
  double* r = box.getAtomPosition(iAtom);
  std::copy(r, r+3, rOld);
  r[0] += 2*stepSize*(random.nextDouble32()-0.5);
  r[1] += 2*stepSize*(random.nextDouble32()-0.5);
  r[2] += 2*stepSize*(random.nextDouble32()-0.5);
  potentialMaster.updateAtom(iAtom);
  numTrials++;
  return true;
}

double MCMoveDisplacement::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOne(iAtom, uNew);
  //printf("uOld %e   uNew %e\n", uOld, uNew);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  chiSum += chi;
  return chi;
}

void MCMoveDisplacement::acceptNotify() {
  //printf("accepted\n");
  potentialMaster.processAtomU(1);
  double uTmp = 0;
  double* r = box.getAtomPosition(iAtom);
  double rSave[3];
  for (int k=0; k<3; k++) {
    rSave[k] = r[k]; r[k] = rOld[k];
  }
  potentialMaster.updateAtom(iAtom);
  potentialMaster.computeOne(iAtom, uTmp);
  for (int k=0; k<3; k++) {
    r[k] = rSave[k];
  }
  potentialMaster.updateAtom(iAtom);
  potentialMaster.processAtomU(-1);
  numAccepted++;
}

void MCMoveDisplacement::rejectNotify() {
  //printf("rejected\n");
  if (iAtom < 0) return;
  double* r = box.getAtomPosition(iAtom);
  std::copy(rOld, rOld+3, r);
  potentialMaster.updateAtom(iAtom);
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveDisplacement::energyChange() {
  return uNew-uOld;
}
