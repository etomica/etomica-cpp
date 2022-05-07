/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-virial.h"

MCMoveChainVirial::MCMoveChainVirial(SpeciesList& sl, Box& b, PotentialMaster& p, Random& r, double s) : MCMove(b,p,r,1), speciesList(sl), sigma(s) {
  tunable = false;
}

MCMoveChainVirial::~MCMoveChainVirial() {}

bool MCMoveChainVirial::doTrial() {
  int nm = box.getTotalNumMolecules();
  if (nm<=1) {
    fprintf(stderr, "Gotta give me more than 1 molecule!\n");
    abort();
  }
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(0, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double rPrev[3];
  double* com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
  std::copy(com, com+3, rPrev);
  for (int iMolecule=1; iMolecule<nm; iMolecule++) {
    double dr[3];
    random.inSphere(dr);
    box.getMoleculeInfo(1, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
    com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
    //printf("move %d %e %e\n", iAtom, r[0]*r[0]+r[1]*r[1]+r[2]*r[2], sigma);
    for (int iAtom=firstAtom; iAtom<=lastAtom; iAtom++) {
      double* r = box.getAtomPosition(iAtom);
      r[0] += (rPrev[0] - com[0]) + dr[0]*sigma;
      r[1] += (rPrev[1] - com[1]) + dr[1]*sigma;
      r[2] += (rPrev[2] - com[2]) + dr[2]*sigma;
    }
  }
  numTrials++;
  numAccepted++;
  return true;
}
