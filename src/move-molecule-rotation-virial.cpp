/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-virial.h"
#include "alloc2d.h"

MCMoveMoleculeRotateVirial::MCMoveMoleculeRotateVirial(SpeciesList& sl, int is, Box& b, PotentialMaster& p, Random& r, double ss, Cluster& c) : MCMove(b,p,r,ss), speciesList(sl), cluster(c), iSpecies(is) {
  Species *sp = speciesList.get(iSpecies);
  int na = sp->getNumAtoms();
  rOld = (double**)malloc2D(na, 3, sizeof(double));
  maxStepSize = M_PI;
}

MCMoveMoleculeRotateVirial::~MCMoveMoleculeRotateVirial() {
  free2D((void**)rOld);
}

bool MCMoveMoleculeRotateVirial::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = box.getNumMolecules(iSpecies);
  iMolecule = box.getGlobalMoleculeIndex(iSpecies, random.nextInt(nm));
  wOld = fabs(cluster.getValues()[0]);
  int axis = random.nextInt(3);
  double theta = stepSize*2*(random.nextDouble32()-0.5);
  mat.setSimpleAxisAngle(axis, theta);
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double* center = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
  double* rfirst = box.getAtomPosition(firstAtom);
  double* rlast = box.getAtomPosition(lastAtom);
  double bvec[] = {rlast[0]-rfirst[0], rlast[1]-rfirst[1], rlast[2]-rfirst[2]};
  for (int iAtom = firstAtom; iAtom <= lastAtom; iAtom++) {
    double* r = box.getAtomPosition(iAtom);
    std::copy(r, r+3, rOld[iAtom-firstAtom]);
    mat.transformAbout(r, center, box);
  }
  bvec[0] = rlast[0]-rfirst[0]; bvec[1] = rlast[1]-rfirst[1]; bvec[2] = rlast[2]-rfirst[2];
  cluster.trialNotify();
  numTrials++;
  return true;
}

double MCMoveMoleculeRotateVirial::getChi(double T) {
  wNew = fabs(cluster.getValues()[0]);
  double chi = wNew>wOld ? 1 : wNew/wOld;
  chiSum += chi;
  return chi;
}

void MCMoveMoleculeRotateVirial::acceptNotify() {
  //printf("accepted\n");
  numAccepted++;
}

void MCMoveMoleculeRotateVirial::rejectNotify() {
  //printf("rejected\n");
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  for (int iAtom = firstAtom; iAtom <= lastAtom; iAtom++) {
    double* r = box.getAtomPosition(iAtom);
    std::copy(rOld[iAtom-firstAtom], rOld[iAtom-firstAtom]+3, r);
  }
}

double MCMoveMoleculeRotateVirial::energyChange() {
  //bogus
  return 0;
}
