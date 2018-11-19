/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-volume.h"

MCMoveVolume::MCMoveVolume(Box& b, PotentialMaster& pm, Random& r, double p, double ss, SpeciesList& sl, Meter& oldPE) : MCMove(b,pm,r,ss), pressure(p), speciesList(sl), oldMeterPE(oldPE) {
  callbacks.push_back(&pce);
}

void MCMoveVolume::scaleVolume(double s) {
  const double* bs = box.getBoxSize();
  box.scaleBoxTo(bs[0]*s, bs[1]*s, bs[2]*s);
}

bool MCMoveVolume::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  uOld = oldMeterPE.getData()[0];
  const double* bs = box.getBoxSize();
  vOld = bs[0]*bs[1]*bs[2];
  lnScale = (2*random.nextDouble32()-1)*stepSize;
  scale = exp(lnScale);
  scaleVolume(scale);
  potentialMaster.updateVolume();
  numTrials++;
  return true;
}

double MCMoveVolume::getChi(double T) {
  pce.reset();
  potentialMaster.computeAll(callbacks);
  uNew = pce.getData()[0];
  double vNew = vOld*scale*scale*scale;
  double deltaH = pressure * (vNew - vOld) + uNew - uOld;
  double x = box.getTotalNumMolecules() * 3*lnScale - deltaH/T;
  double chi = x>0 ? 1 : exp(x);
  chiSum += chi;
  return chi;
}

void MCMoveVolume::acceptNotify() {
  //printf("accepted => %f\n", box.getBoxSize()[0]);
  numAccepted++;
}

void MCMoveVolume::rejectNotify() {
  scaleVolume(1/scale);
  //printf("rejected => %f\n", box.getBoxSize()[0]);
  potentialMaster.updateVolume();
  // we have to recompute energy so that atom energies are re-updated
  pce.reset();
  potentialMaster.computeAll(callbacks);
}

double MCMoveVolume::energyChange() {
  return uNew-uOld;
}
