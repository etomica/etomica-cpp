/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "cluster.h"

ClusterFlipped::ClusterFlipped(Cluster& c, SpeciesList& sl, Box& b, bool cached) : Cluster(b.getTotalNumMolecules(),c.numValues(),cached), wrappedCluster(c), speciesList(sl), box(b) {
  flippedAtoms = (bool*)malloc(numMolecules*sizeof(bool));
}

ClusterFlipped::~ClusterFlipped() {
  free(flippedAtoms);
}

void ClusterFlipped::flip(int iMolecule) {
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double* com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
  for (int iAtom=firstAtom; iAtom<=lastAtom; iAtom++) {
    double* r = box.getAtomPosition(iAtom);
    for (int i=0; i<3; i++) {
      r[i] = 2*com[i] - r[i];
    }
  }
}

const double* ClusterFlipped::getValues() {
  if (useCache && !cacheDirty) {
    return values;
  }
  cacheDirty = false;

  bool flipit = false;

  double minFlipDistance = 20;
  double minR2 = minFlipDistance*minFlipDistance;
  for(int iMol1=0; iMol1<numMolecules; iMol1++){
    flippedAtoms[iMol1] = false;
    int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
    box.getMoleculeInfo(iMol1, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
    double* com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
    double com1[3] = {com[0],com[1],com[2]};
    for (int iMol2=iMol1+1; iMol2<numMolecules; iMol2++) {
      box.getMoleculeInfo(iMol2, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
      com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
      double dr[3] = {com[0]-com1[0], com[1]-com1[1], com[2]-com1[2]};
      if (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] > minR2) {
        flipit = true;
      }
    }
  }

  const double* foo = wrappedCluster.getValues();
  for (int i=0; i<nValues; i++) {
    values[i] = foo[i];
  }
  if (!flipit) {
    return values;
  }

  while (true) {
    bool didFlipTrue = false;
    for (int i=0; !didFlipTrue && i<numMolecules; i++) {
      flippedAtoms[i] = !flippedAtoms[i];
      didFlipTrue = flippedAtoms[i];
      flip(i);
    }
    if (!didFlipTrue) {
      // if we flipped every molecule from true to false, we must be done
      break;
    }
    const double* foo = wrappedCluster.getValues();
    for (int i=0; i<nValues; i++) {
      values[i] += foo[i];
    }
  }

  for  (int i=0; i<nValues; i++) {
    values[i] /= pow(2, numMolecules);
  }

  return values;
}
