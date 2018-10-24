/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move.h"
#include "alloc2d.h"

MCMoveMoleculeDisplacement::MCMoveMoleculeDisplacement(Box& b, PotentialMaster& p, Random& r, double ss) : MCMove(b,p,r,ss), numOldPositions(0), oldPositions(nullptr), iSpecies(-1) {
}

MCMoveMoleculeDisplacement::~MCMoveMoleculeDisplacement() {
  free2D((void**)oldPositions);
}

bool MCMoveMoleculeDisplacement::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = iSpecies<0 ? box.getTotalNumMolecules() : box.getNumMolecules(iSpecies);
  if (nm==0) {
    iMolecule = -1;
    return false;
  }
  iMolecule = random.nextInt(nm);
  int iMoleculeInSpecies;
  if (iSpecies>=0) {
    iMoleculeInSpecies = iMolecule;
    iMolecule = box.getGlobalMoleculeIndex(iSpecies, iMolecule);
  }
  uOld = potentialMaster.oldMoleculeEnergy(iMolecule);
  deltaR[0] = 2*stepSize*(random.nextDouble32()-0.5);
  deltaR[1] = 2*stepSize*(random.nextDouble32()-0.5);
  deltaR[2] = 2*stepSize*(random.nextDouble32()-0.5);
  int na = numOldPositions;
  if (iSpecies < 0) {
    int mySpecies;
    box.getMoleculeInfo(iMolecule, mySpecies, iMoleculeInSpecies, iAtomFirst, iAtomLast);
    na = iAtomLast-iAtomFirst+1;
    if (na>numOldPositions) {
      oldPositions = (double**)realloc2D((void**)oldPositions, na, 3, sizeof(double));
      numOldPositions = na;
    }
  }
  else {
    iAtomFirst = box.getFirstAtom(iSpecies, iMoleculeInSpecies);
    iAtomLast = iAtomFirst + na - 1;
  }
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) {
      oldPositions[i][j] = ri[j];
      ri[j] += deltaR[j];
    }
    box.nearestImage(ri);
    potentialMaster.updateAtom(iAtom);
  }
  numTrials++;
  return true;
}

double MCMoveMoleculeDisplacement::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOneMolecule(iMolecule, uNew);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  chiSum += chi;
  return chi;
}

void MCMoveMoleculeDisplacement::acceptNotify() {
  //printf("accepted\n");
  potentialMaster.processAtomU(1);
  int na = iAtomLast-iAtomFirst+1;
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) ri[j] = oldPositions[i][j];
    potentialMaster.updateAtom(iAtom);
  }
  double uTmp = 0;
  potentialMaster.computeOneMolecule(iMolecule, uTmp);
  potentialMaster.processAtomU(-1);
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) ri[j] += deltaR[j];
    box.nearestImage(ri);
    potentialMaster.updateAtom(iAtom);
  }
  numAccepted++;
}

void MCMoveMoleculeDisplacement::rejectNotify() {
  //printf("rejected\n");
  if (iMolecule < 0) return;
  int na = iAtomLast-iAtomFirst+1;
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) ri[j] = oldPositions[i][j];
    potentialMaster.updateAtom(iAtom);
  }
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveMoleculeDisplacement::energyChange() {
  return uNew-uOld;
}

void MCMoveMoleculeDisplacement::setSpecies(int s) {
  iSpecies = s;
  if (s<0) return;
  int numAtoms = box.getSpeciesList().get(s)->getNumAtoms();
  oldPositions = (double**)realloc2D((void**)oldPositions, numAtoms, 3, sizeof(double));
  numOldPositions = numAtoms;
}
