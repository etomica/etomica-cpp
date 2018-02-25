#include "move.h"

MCMoveInsertDelete::MCMoveInsertDelete(Box& b, PotentialMaster& p, Random& r, double m, int s) : MCMove(b,p,r,0), mu(m), iSpecies(s) {
  tunable = false;
}

MCMoveInsertDelete::~MCMoveInsertDelete() {}

bool MCMoveInsertDelete::doTrial() {
  doInsert = random.nextInt(2) == 0;
  int n = box.getNumAtoms();
  if (doInsert) {
    uOld = 0;
    double *bs = box.getBoxSize();
    for (int j=0; j<3; j++) {
      rNew[j] = bs[j]*(random.nextDouble32()-0.5);
    }
    potentialMaster.computeOne(n, rNew, uNew, true);
  }
  else {
    // delete
    if (n>0) {
      iAtom = random.nextInt(n);
      uOld = potentialMaster.oldEnergy(iAtom);
      uNew = 0;
    }
  }
  numTrials++;
  return true;
}

double MCMoveInsertDelete::getChi(double T) {
  int n = box.getNumAtoms();
  if (n==0 && !doInsert) return 0;

  double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];

  double x = -(uNew-uOld);
  double a;
  if (doInsert) {
    x += mu;
    a = vol/(n+1);
  }
  else {
    x -= mu;
    a = (n-1)/vol;
  }
  double chi = a*exp(x/T);
  //printf("%d %d %f %f %f %f %f %f\n", doInsert, n, a, uOld, uNew, mu, x, chi);
  if (chi>1) chi = 1;
  chiSum += chi;
  return chi;
}

void MCMoveInsertDelete::acceptNotify() {
  int n = box.getTotalNumMolecules();
  if (doInsert) {
    //printf("accept insert %d\n", n+1);
    box.setNumMolecules(iSpecies, n+1);
    double *ri = box.getAtomPosition(n);
    std::copy(rNew, rNew+3, ri);
    potentialMaster.newAtom();
    potentialMaster.processAtomU(1);
  }
  else {
    //printf("accept delete %d\n", n-1);
    double *ri = box.getAtomPosition(iAtom);
    double uTmp = 0;
    potentialMaster.computeOne(iAtom, ri, uTmp, false);
    potentialMaster.processAtomU(-1);
    // this removes iAtom (updates cell lists) and then
    // moves the last atom into its position
    potentialMaster.removeAtom(iAtom);
    double *rn = box.getAtomPosition(n-1);
    std::copy(rn, rn+3, ri);
    box.setNumMolecules(iSpecies, n-1);
  }

  numAccepted++;
}

void MCMoveInsertDelete::rejectNotify() {
  //if (doInsert) printf("reject insert\n");
  //else printf("reject delete\n");
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveInsertDelete::energyChange() {
  return uNew-uOld;
}
