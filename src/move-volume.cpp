/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-volume.h"

MCMoveVolume::MCMoveVolume(Box& b, PotentialMaster& pm, Random& r, double p, double ss, SpeciesList& sl, Meter& oldPE) : MCMove(b,pm,r,ss), pressure(p), speciesList(sl), oldMeterPE(oldPE) {
  callbacks.push_back(&pce);
}

void MCMoveVolume::scaleVolume(double s) {
  int nm = box.getTotalNumMolecules();
  s -= 1;
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* iPos = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    double dr[3] = {iPos[0]*s, iPos[1]*s, iPos[2]*s};
    for (int jAtom=iFirstAtom; jAtom<=iLastAtom; jAtom++) {
      double* rj = box.getAtomPosition(jAtom);
      for (int k=0; k<3; k++) rj[k] += dr[k];
    }
  }
  const double* bs = box.getBoxSize();
  s += 1;
  double bx = bs[0]*s, by = bs[1]*s, bz = bs[2]*s;
  box.setBoxSize(bx, by, bz);
  int na = box.getNumAtoms();
  for (int iAtom=0; iAtom<na; iAtom++) {
    double* ri = box.getAtomPosition(iAtom);
    box.nearestImage(ri);
  }
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
