#include "move.h"
#include "alloc2d.h"

MCMoveMoleculeDisplacement::MCMoveMoleculeDisplacement(Box& b, PotentialMaster& p, Random& r, double ss) : MCMove(b,p,r,ss), numOldPositions(0), oldPositions(nullptr) {
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
  int na = iAtomLast-iAtomFirst+1;
  if (na>numOldPositions) {
    oldPositions = (double**)realloc2D((void**)oldPositions, na, 3, sizeof(double));
    numOldPositions = na;
  }
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);

    for (int j=0; j<3; j++) {
      oldPositions[i][j] = ri[j];
      ri[j] += deltaR[j];
    }
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
  int na = iAtomLast-iAtomFirst+1;
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) ri[j] = oldPositions[i][j];
  }
  double uTmp = 0;
  potentialMaster.computeOneMolecule(iMolecule, uTmp, false);
  potentialMaster.processAtomU(-1);
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) ri[j] += deltaR[j];
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
  }
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveMoleculeDisplacement::energyChange() {
  return uNew-uOld;
}
