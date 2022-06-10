/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "meter.h"
#include "potential-master.h"
#include "matrix.h"

MeterWidomInsertion::MeterWidomInsertion(Box& b, int is, PotentialMaster& pm, Random& rand, double t, int n) : Meter(1), box(b), iSpecies(is), potentialMaster(pm), random(rand), temperature(t), numTrials(n) {
}

MeterWidomInsertion::~MeterWidomInsertion() {}

double* MeterWidomInsertion::getData() {
  Matrix* h = box.getH();
  int N = box.getNumMolecules(iSpecies);
  box.setNumMolecules(iSpecies, N+1);
  potentialMaster.newMolecule(iSpecies);
  int iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(N, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);

  double sum = 0;
  for (int iTrial=0; iTrial<numTrials; iTrial++) {
    double *r = box.getAtomPosition(firstAtom); // let's just pretend that's OK
    for (int k=0; k<3; k++) r[k] = random.nextDouble() - 0.5;
    h->transform(r);
    for (int iAtom=firstAtom; iAtom<=lastAtom; iAtom++) potentialMaster.updateAtom(iAtom);
    double u = 0;
    potentialMaster.computeOneMolecule(N, u);
    sum += exp(-u/temperature);
  }
  potentialMaster.unsetAtomDU();
  box.setNumMolecules(iSpecies, N);
  data[0] = sum/numTrials;
  return data;
}
