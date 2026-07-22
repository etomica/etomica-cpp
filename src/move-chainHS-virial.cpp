/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "move-virial.h"

MCMoveChainVirial::MCMoveChainVirial(SpeciesList& sl, Box& b, Random& r, double s) : MCMove(b,r,1), speciesList(sl), sigma(s) {
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
  if (!fixedCOM) {
    std::copy(box.getAtomPosition(firstAtom), box.getAtomPosition(firstAtom) + 3, com);
  }

  // double* C0 = box.getAtomPosition(firstAtom);
  // printf("molecule 0 C %f %f %f ", C0[0], C0[1], C0[2]);

  // if (fabs(com[2]) > 0.1) exit(0);
  std::copy(com, com+3, rPrev);
  for (int iMolecule=1; iMolecule<nm; iMolecule++) {
    double dr[3];
    // random.inSphere(dr);
    double discreteStepSize = 0.05;
    long n = sigma / discreteStepSize;
    long pSum = n * (n+1) * (1+2*n) / 6;
    long iSum = 0;
    int i;
    double rand = random.nextDouble();
    for (i=n;i>0;i--) {
      iSum += i*i;
      if (iSum>=pSum*rand) break;
    }
    dr[0] = i*discreteStepSize/sigma;
    // printf("%f\n", dr[0]*sigma);
    // dr[0] = random.nextDouble()*2-1;
    dr[1] = dr[2] = 0;
    double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
    com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
    if (!fixedCOM) {
      std::copy(box.getAtomPosition(firstAtom), box.getAtomPosition(firstAtom) + 3, com);
    }

    for (int k=0; k<3; k++) {
      dr[k] = -com[k] + rPrev[k] + dr[k]*sigma;

    }
    // printf("move %d %e %e\n", iMolecule, dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2], sigma);
    // printf("%f %f %f \n", com[0], com[1], com[2]);

    for (int iAtom=firstAtom; iAtom<=lastAtom; iAtom++) {
      double* r = box.getAtomPosition(iAtom);

      r[0] += dr[0];
      r[1] += dr[1];
      r[2] += dr[2];
      // if (iAtom==firstAtom)
      //  {
      //    printf(" molecule 1 C i: %d %f: %f %f %f\n", i, rand, r[0], r[1], r[2]);
      //  }
// }

    }
    for (int k=0; k<3; k++) {
      rPrev[k] = com[k] + dr[k];
    }
    // com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
    // printf("%f %f %f \n", com[0], com[1], com[2]);

  }

  numTrials++;
  numAccepted++;
  return true;
}
