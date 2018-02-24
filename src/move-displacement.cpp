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
  box.nearestImage(r);
  numTrials++;
  return true;
}

double MCMoveDisplacement::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOne(iAtom, box.getAtomPosition(iAtom), uNew, true);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  chiSum += chi;
  return chi;
}

void MCMoveDisplacement::acceptNotify() {
  //printf("accepted\n");
  potentialMaster.processAtomU(1);
  double uTmp = 0;
  potentialMaster.computeOne(iAtom, rOld, uTmp, false);
  potentialMaster.processAtomU(-1);
  potentialMaster.updateAtom(iAtom);
  numAccepted++;
}

void MCMoveDisplacement::rejectNotify() {
  //printf("rejected\n");
  if (iAtom < 0) return;
  double* r = box.getAtomPosition(iAtom);
  std::copy(rOld, rOld+3, r);
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveDisplacement::energyChange() {
  return uNew-uOld;
}
