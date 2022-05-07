/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-virial.h"
#include "alloc2d.h"

MCMoveMoleculeDisplacementVirial::MCMoveMoleculeDisplacementVirial(SpeciesList& sl, int is, Box& b, PotentialMaster& p, Random& r, double ss, Cluster& c) : MCMove(b,p,r,ss), cluster(c), iSpecies(is) {
  Species *sp = sl.get(iSpecies);
  int na = sp->getNumAtoms();
  if (na<=1 && iSpecies==0) {
    fprintf(stderr, "Cannot translate molecule 0");
    exit(1);
  }
  rOld = (double**)malloc2D(na, 3, sizeof(double));
}

MCMoveMoleculeDisplacementVirial::~MCMoveMoleculeDisplacementVirial() {
  free2D((void**)rOld);
}

bool MCMoveMoleculeDisplacementVirial::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = box.getNumMolecules(iSpecies);
  iMolecule = 0;
  while (iMolecule == 0) {
    iMolecule = box.getGlobalMoleculeIndex(iSpecies, random.nextInt(nm));
  }
  wOld = fabs(cluster.getValues()[0]);
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double deltaR[3];
  deltaR[0] = 2*stepSize*(random.nextDouble32()-0.5);
  deltaR[1] = 2*stepSize*(random.nextDouble32()-0.5);
  deltaR[2] = 2*stepSize*(random.nextDouble32()-0.5);
  for (int iAtom = firstAtom; iAtom <= lastAtom; iAtom++) {
    double* r = box.getAtomPosition(iAtom);
    std::copy(r, r+3, rOld[iAtom-firstAtom]);
    r[0] += deltaR[0];
    r[1] += deltaR[1];
    r[2] += deltaR[2];
  }
  cluster.trialNotify();
  numTrials++;
  return true;
}

double MCMoveMoleculeDisplacementVirial::getChi(double T) {
  wNew = fabs(cluster.getValues()[0]);
  double chi = wNew>wOld ? 1 : wNew/wOld;
  chiSum += chi;
  return chi;
}

void MCMoveMoleculeDisplacementVirial::acceptNotify() {
  //printf("accepted\n");
  numAccepted++;
}

void MCMoveMoleculeDisplacementVirial::rejectNotify() {
  //printf("rejected\n");
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  for (int iAtom = firstAtom; iAtom <= lastAtom; iAtom++) {
    double* r = box.getAtomPosition(iAtom);
    std::copy(rOld[iAtom-firstAtom], rOld[iAtom-firstAtom]+3, r);
  }
}

double MCMoveMoleculeDisplacementVirial::energyChange() {
  //bogus
  return 0;
}
