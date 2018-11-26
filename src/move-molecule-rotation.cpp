/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move.h"
#include "alloc2d.h"

MCMoveMoleculeRotate::MCMoveMoleculeRotate(SpeciesList& sl, Box& b, PotentialMaster& p, Random& r) : MCMove(b,p,r,0.5), speciesList(sl), numOldPositions(0), oldPositions(nullptr), iSpecies(-1) {
  maxStepSize = M_PI;
}

MCMoveMoleculeRotate::~MCMoveMoleculeRotate() {
  free2D((void**)oldPositions);
}

void MCMoveMoleculeRotate::transformAbout(double* r, const double* rOld, double* center) {
  mat.transformAbout(r, center, box);
  // ugh.  that transformed and wrapped into the box.
  // we need to (maybe) unwrap to be on the same side of the box it started at
  double dr[3];
  for (int k=0; k<3; k++) dr[k] = r[k] - rOld[k];
  box.nearestImage(dr);
  for (int k=0; k<3; k++) r[k] = rOld[k] + dr[k];
}

bool MCMoveMoleculeRotate::doTrial() {
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
  if (false) {
    // see if we can recompute the molecule energy and get the same result
    // this messes up future computation, so abort afterwards
    printf("got old, now recompute\n");
    double uTmp = 0;
    potentialMaster.computeOneMolecule(iMolecule, uTmp);
    printf("%d uOlds %f %f\n", iMolecule, uTmp, uOld);
    abort();
  }
  int na = numOldPositions;
  mySpecies = iSpecies;
  if (iSpecies<0) {
    box.getMoleculeInfo(iMolecule, mySpecies, iMoleculeInSpecies, iAtomFirst, iAtomLast);
    na = iAtomLast-iAtomFirst+1;
    if (na==1) {
      iMolecule = -1;
      return false;
    }
    if (na>numOldPositions) {
      oldPositions = (double**)realloc2D((void**)oldPositions, na, 3, sizeof(double));
      numOldPositions = na;
    }
  }
  else {
    iAtomFirst = box.getFirstAtom(mySpecies, iMoleculeInSpecies);
    iAtomLast = iAtomFirst + na - 1;
  }
  int axis = random.nextInt(3);
  double theta = stepSize*2*(random.nextDouble32()-0.5);
  mat.setSimpleAxisAngle(axis, theta);
  Species* species = speciesList.get(mySpecies);
  double* center = species->getMoleculeCOM(box, iAtomFirst, iAtomLast);
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    for (int j=0; j<3; j++) {
      oldPositions[i][j] = ri[j];
    }
    transformAbout(ri, oldPositions[i], center);
    // and then rewrap into the box if the potentialMaster wants it
    potentialMaster.updateAtom(iAtom);
  }
  numTrials++;
  return true;
}

double MCMoveMoleculeRotate::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOneMolecule(iMolecule, uNew);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  //printf("%f => %f  ==>  %f\n", uOld, uNew, uNew-uOld);
  chiSum += chi;
  return chi;
}

void MCMoveMoleculeRotate::acceptNotify() {
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
  // this call is designed to set up the next call.  uTmp won't necessarily be correct
  potentialMaster.processAtomU(-1);
  Species* species = speciesList.get(mySpecies);
  double* center = species->getMoleculeCOM(box, iAtomFirst, iAtomLast);
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);
    transformAbout(ri, oldPositions[i], center);
    potentialMaster.updateAtom(iAtom);
  }
  numAccepted++;
}

void MCMoveMoleculeRotate::rejectNotify() {
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

double MCMoveMoleculeRotate::energyChange() {
  return uNew-uOld;
}

void MCMoveMoleculeRotate::setSpecies(int s) {
  iSpecies = s;
  if (s<0) return;
  int numAtoms = box.getSpeciesList().get(s)->getNumAtoms();
  if (numAtoms < 2) {
    fprintf(stderr, "Need at least 2 atoms in a molecule for rotation!\n");
    abort();
  }
  oldPositions = (double**)realloc2D((void**)oldPositions, numAtoms, 3, sizeof(double));
  numOldPositions = numAtoms;
}
