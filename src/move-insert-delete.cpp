#include "move.h"

MCMoveInsertDelete::MCMoveInsertDelete(Box& b, PotentialMaster& p, Random& r, double m, int s) : MCMove(b,p,r,0), mu(m), iSpecies(s), numAtoms(box.getSpeciesList().get(iSpecies)->getNumAtoms()) {
  tunable = false;
}

MCMoveInsertDelete::~MCMoveInsertDelete() {}

bool MCMoveInsertDelete::doTrial() {
  doInsert = random.nextInt(2) == 0;
  int n = box.getNumMolecules(iSpecies);
  if (doInsert) {
    xMolecule = n;
    box.setNumMolecules(iSpecies, n+1);
    uOld = 0;
    const double *bs = box.getBoxSize();
    double mPos[3];
    for (int k=0; k<3; k++) {
      mPos[k] = bs[k]*(random.nextDouble32()-0.5);
    }
    firstAtom = box.getFirstAtom(iSpecies, n);
    if (numAtoms==1) {
      double *rj = box.getAtomPosition(firstAtom);
      for (int k=0; k<3; k++) rj[k] = mPos[k];
      iMolecule = box.getGlobalMoleculeIndex(iSpecies, n);
      potentialMaster.newMolecule(iSpecies);
      potentialMaster.computeOne(firstAtom, uNew);
    }
    else {
      double *r0;
      rotMat.randomize(random);
      for (int j=0; j<numAtoms; j++) {
        int jAtom = firstAtom+j;
        double *rj = box.getAtomPosition(jAtom);
        for  (int k=0; k<3; k++) rj[k] += mPos[k];
        if (j==0) r0 = rj;
        else {
          rotMat.transformAbout(rj, r0, box);
        }
      }
      iMolecule = box.getGlobalMoleculeIndex(iSpecies, n);
      potentialMaster.newMolecule(iSpecies);
      potentialMaster.computeOneMolecule(iMolecule, uNew);
    }
  }
  else {
    // delete
    if (n>0) {
      xMolecule = random.nextInt(box.getNumMolecules(iSpecies));
      iMolecule = box.getGlobalMoleculeIndex(iSpecies, xMolecule);
      uOld = potentialMaster.oldMoleculeEnergy(iMolecule);
      uNew = 0;
    }
    else {
      numTrials++;
      return false;
    }
  }
  numTrials++;
  return true;
}

double MCMoveInsertDelete::getChi(double T) {
  int n = box.getNumMolecules(iSpecies);

  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];

  double x = -(uNew-uOld);
  double a;
  if (doInsert) {
    x += mu;
    a = vol/n;
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
  if (doInsert) {
    //printf("accept insert %d\n", iMolecule);
    potentialMaster.processAtomU(1);
  }
  else {
    //printf("accept delete %d\n", iMolecule);
    if (numAtoms==1) {
      double uTmp;
      firstAtom = box.getFirstAtom(iSpecies, xMolecule);
      potentialMaster.computeOne(firstAtom, uTmp);
      potentialMaster.processAtomU(-1);
    }
    else {
      double uTmp;
      potentialMaster.computeOneMolecule(iMolecule, uTmp);
      potentialMaster.processAtomU(-1);
    }
    // this removes iMolecule (updates cell lists) and then
    // moves the last atom into its position
    potentialMaster.removeMolecule(iSpecies, xMolecule);

    int jMolecule = box.getNumMolecules(iSpecies)-1;
    if (jMolecule>xMolecule) {
      int xFirstAtom = box.getFirstAtom(iSpecies, xMolecule);
      int jFirstAtom = box.getFirstAtom(iSpecies, jMolecule);
      for (int k=0; k<numAtoms; k++) {
        int kxAtom = xFirstAtom + k;
        int kjAtom = jFirstAtom + k;
        double *rx = box.getAtomPosition(kxAtom);
        double *rj = box.getAtomPosition(kjAtom);
        for (int l=0; l<3; l++) rx[l] = rj[l];
      }
    }
    box.setNumMolecules(iSpecies, jMolecule);
  }

  numAccepted++;
}

void MCMoveInsertDelete::rejectNotify() {
  //if (doInsert) printf("reject insert\n");
  //else printf("reject delete\n");
  if (doInsert) {
    potentialMaster.removeMolecule(iSpecies, iMolecule);
    box.setNumMolecules(iSpecies, xMolecule);
  }
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveInsertDelete::energyChange() {
  return uNew-uOld;
}
