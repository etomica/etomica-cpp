/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-virial.h"
#include "alloc2d.h"

MCMoveMoleculeDisplacementVirial::MCMoveMoleculeDisplacementVirial(SpeciesList& sl, int is, Box& b, Random& r, double ss, Cluster& c, double dc, double dss) :
  MCMove(b, r, ss), cluster(c), iSpecies(is), discreteCutOff(dc), discreteStepSize(dss)
{
  Species* sp = sl.get(iSpecies);
  int na = sp->getNumAtoms();
  rOld = (double**)malloc2D(na, 3, sizeof(double));
  for (int i = 0; i <= 90; i++)
  {
    pisum[i] = hcount[i] = 0;
  }
  maxStepSize = 100;
}

MCMoveMoleculeDisplacementVirial::~MCMoveMoleculeDisplacementVirial() {
  free2D((void**)rOld);
}

bool MCMoveMoleculeDisplacementVirial::doTrial() {
  if (false) return false;
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = box.getNumMolecules(iSpecies);
  iMolecule = 0;
  while (iMolecule == 0) {
    iMolecule = box.getGlobalMoleculeIndex(iSpecies, random.nextInt(nm));
  }
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double rTemp = discreteStepSize > 0 ? box.getAtomPosition(firstAtom)[0] : 1;
  wOld = fabs(cluster.getValues()[0])*rTemp*rTemp;

  double deltaR[3];
  if (discreteStepSize > 0)
  {
    int n = stepSize /  discreteStepSize;
    int iStep = random.nextInt(2*n) - n;
    if (iStep >  -1) iStep++;
    deltaR[0] = discreteStepSize*iStep;
    deltaR[1] = 0*stepSize*(random.nextDouble32()-0.5);
    deltaR[2] = 0*stepSize*(random.nextDouble32()-0.5);

  }
  else
  {
    deltaR[0] = stepSize*(random.nextDouble32()-0.5);;
    deltaR[1] = stepSize*(random.nextDouble32()-0.5);
    deltaR[2] = stepSize*(random.nextDouble32()-0.5);
    // printf("%f %f \n", stepSize, deltaR[0]);

  }
  for (int iAtom = firstAtom; iAtom <= lastAtom; iAtom++) {
    double* r = box.getAtomPosition(iAtom);
    std::copy(r, r+3, rOld[iAtom-firstAtom]);
    r[0] += deltaR[0];
    r[1] += deltaR[1];
    r[2] += deltaR[2];
    // if (iAtom == firstAtom) printf("%f \n", r[0]);
    // if (iAtom == firstAtom && r[0] > 3)
    // {
    //   printf("breakpoint \n");
    // }
  }
  rNew = discreteStepSize > 0 ? box.getAtomPosition(firstAtom)[0] : 1;
  cluster.trialNotify();
  numTrials++;
  return true;
}

double MCMoveMoleculeDisplacementVirial::getChi(double T) {
  wNew = fabs(cluster.getValues()[0])*rNew*rNew;
  double chi = wNew>wOld ? 1 : wNew/wOld;
  // printf("%f %f %f %f\n", wNew, wOld, chi, discreteCutOff);
  if (rNew > discreteCutOff) chi = 0;
  chiSum += chi;
  return chi;
}

void MCMoveMoleculeDisplacementVirial::acceptNotify() {
  // printf("accepted\n");
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
  // printf("rejected\n");
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
