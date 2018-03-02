#include "move.h"

MCMoveMoleculeDisplacement::MCMoveMoleculeDisplacement(Box& b, PotentialMaster& p, Random& r, double ss) : MCMove(b,p,r,ss) {
}

MCMoveMoleculeDisplacement::~MCMoveMoleculeDisplacement() {}

bool MCMoveMoleculeDisplacement::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = box.getTotalNumMolecules();
  if (nm==0) {
    iMolecule = -1;
    return false;
  }
  iMolecule = random.nextInt(nm);
  uOld = potentialMaster.oldMoleculeEnergy(iMolecule);
  deltaR[0] = 2*stepSize*(random.nextDouble32()-0.5);
  deltaR[1] = 2*stepSize*(random.nextDouble32()-0.5);
  deltaR[2] = 2*stepSize*(random.nextDouble32()-0.5);
  int iSpecies;
  box.getMoleculeInfo(iMolecule, iSpecies, iAtomFirst, iAtomLast);
  for (int iAtom=iAtomFirst; iAtom<=iAtomLast; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) ri[j] += deltaR[j];
    box.nearestImage(ri);
  }
  numTrials++;
  return true;
}

double MCMoveMoleculeDisplacement::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOneMolecule(iMolecule, uNew, true);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  chiSum += chi;
  return chi;
}

void MCMoveMoleculeDisplacement::acceptNotify() {
  //printf("accepted\n");
  potentialMaster.processAtomU(1);
  double uTmp = 0;
  potentialMaster.computeOneMolecule(iMolecule, uTmp, false);
  potentialMaster.processAtomU(-1);
  for (int iAtom=iAtomFirst; iAtom<=iAtomLast; iAtom++) {
    potentialMaster.updateAtom(iAtom);
  }
  numAccepted++;
}

void MCMoveMoleculeDisplacement::rejectNotify() {
  //printf("rejected\n");
  if (iMolecule < 0) return;
  for (int iAtom=iAtomFirst; iAtom<=iAtomLast; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) ri[j] -= deltaR[j];
    box.nearestImage(ri);
  }
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveMoleculeDisplacement::energyChange() {
  return uNew-uOld;
}
