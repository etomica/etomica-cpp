#include "move.h"
#include "alloc2d.h"

MCMoveMoleculeRotate::MCMoveMoleculeRotate(AtomInfo& ai, Box& b, PotentialMaster& p, Random& r) : MCMove(b,p,r,0.5), atomInfo(ai), numOldPositions(0), oldPositions(nullptr) {
}

MCMoveMoleculeRotate::~MCMoveMoleculeRotate() {
  free2D((void**)oldPositions);
}

bool MCMoveMoleculeRotate::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = box.getTotalNumMolecules();
  if (nm==0) {
    iMolecule = -1;
    return false;
  }
  iMolecule = random.nextInt(nm);
  int iSpecies, iMoleculeInSpecies;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iAtomFirst, iAtomLast);
  if (iAtomFirst==iAtomLast) {
    iMolecule = -1;
    return false;
  }
  uOld = potentialMaster.oldMoleculeEnergy(iMolecule);
  int na = iAtomLast-iAtomFirst+1;
  if (na>numOldPositions) {
    oldPositions = (double**)realloc2D((void**)oldPositions, na, 3, sizeof(double));
    numOldPositions = na;
  }
  int axis = random.nextInt(3);
  double theta = stepSize*random.nextDouble32();
  mat.setSimpleAxisAngle(axis, theta);
  double center[3] = {0,0,0};
  double totMass = 0;
  double *r0 = box.getAtomPosition(iAtomFirst);
  for (int i=0; i<na; i++) {
    int iAtom=iAtomFirst+i;
    double *ri = box.getAtomPosition(iAtom);
    int iType = box.getAtomType(iAtom);
    double mass = atomInfo.getMass(iType);
    totMass += mass;
    if (i==0) {
      for (int k=0; k<3; k++) center[k] = mass*ri[k];
    }
    else {
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = ri[k] - r0[k];
      box.nearestImage(dr);
      for (int k=0; k<3; k++) center[k] += mass*(r0[k] + dr[k]);
    }
  }
  for (int k=0; k<3; k++) center[k] /= totMass;
  for (int i=0; i<na; i++) {
    int iAtom=iAtomFirst+i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) {
      oldPositions[i][j] = ri[j];
    }
    mat.transformAbout(ri, center, box);
    potentialMaster.updateAtom(iAtom);
  }
  numTrials++;
  return true;
}

double MCMoveMoleculeRotate::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOneMolecule(iMolecule, uNew);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  chi = 0;
  chiSum += chi;
  return chi;
}

void MCMoveMoleculeRotate::acceptNotify() {
  //printf("accepted\n");
  potentialMaster.processAtomU(1);
  int na = iAtomLast-iAtomFirst+1;
  for (int i=0; i<na; i++) {
    int iAtom=iAtomFirst+i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) {
      ri[j] = oldPositions[i][j];
    }
    potentialMaster.updateAtom(iAtom);
  }
  double uTmp = 0;
  potentialMaster.computeOneMolecule(iMolecule, uTmp);
  potentialMaster.processAtomU(-1);
  double *r0 = box.getAtomPosition(iAtomFirst);
  for (int i=0; i<na; i++) {
    int iAtom=iAtomFirst+i;
    double *ri = box.getAtomPosition(iAtom);
    mat.transformAbout(ri, r0, box);
    potentialMaster.updateAtom(iAtom);
  }
  numAccepted++;
}

void MCMoveMoleculeRotate::rejectNotify() {
  //printf("rejected\n");
  if (iMolecule < 0) return;
  int na = iAtomLast-iAtomFirst+1;
  for (int i=0; i<na; i++) {
    int iAtom=iAtomFirst+i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) {
      ri[j] = oldPositions[i][j];
    }
    potentialMaster.updateAtom(iAtom);
  }
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveMoleculeRotate::energyChange() {
  return uNew-uOld;
}
