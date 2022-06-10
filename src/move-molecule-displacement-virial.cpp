/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-virial.h"
#include "alloc2d.h"

MCMoveMoleculeDisplacementVirial::MCMoveMoleculeDisplacementVirial(SpeciesList& sl, int is, Box& b, PotentialMaster& p, Random& r, double ss, Cluster& c) : MCMove(b,p,r,ss), cluster(c), iSpecies(is) {
  Species *sp = sl.get(iSpecies);
  int na = sp->getNumAtoms();
  rOld = (double**)malloc2D(na, 3, sizeof(double));
  for (int i=0; i<=90; i++) {
    pisum[i] = hcount[i] =0;
  }
  maxStepSize = 100;
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
  addToHistogram(wNew);
}

double* MCMoveMoleculeDisplacementVirial::getHistogramPi() {
  for (int i=0; i<=90; i++) {
    if (hcount[i] == 0) piHist[i] = nan("");
    else piHist[i] = pisum[i] / (double)hcount[i];
  }
  return piHist;
}

double* MCMoveMoleculeDisplacementVirial::getHistogram() {
  long total = 0;
  for (int i=0; i<=90; i++) {
    total += hcount[i];
  }

  for (int i=0; i<=90; i++) {
    histogram[i] = hcount[i] / (double)total * 10;
  }
  return histogram;
}

void MCMoveMoleculeDisplacementVirial::addToHistogram(double pi) {
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(1, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double* ra1 = box.getAtomPosition(firstAtom);
  double r2 = ra1[0]*ra1[0] + ra1[1]*ra1[1] + ra1[2]*ra1[2];
  double r1 = sqrt(r2);
  if (r1>1) r1 = 1+log(r1);
  int ridx = (int)(r1*10);
  if (ridx <= 90) {
    pisum[ridx] += pi;
    hcount[ridx]++;
  }
}

void MCMoveMoleculeDisplacementVirial::rejectNotify() {
  //printf("rejected\n");
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  for (int iAtom = firstAtom; iAtom <= lastAtom; iAtom++) {
    double* r = box.getAtomPosition(iAtom);
    std::copy(rOld[iAtom-firstAtom], rOld[iAtom-firstAtom]+3, r);
  }
  addToHistogram(wOld);
}

double MCMoveMoleculeDisplacementVirial::energyChange() {
  //bogus
  return 0;
}
