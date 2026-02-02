/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "cluster.h"

ClusterFlippedPoints::ClusterFlippedPoints(Cluster& c, SpeciesList& sl, Box& b, bool cached, int mfd) :
  Cluster(b.getTotalNumMolecules(), c.numValues(), cached), wrappedCluster(c), speciesList(sl), box(b),
  minFlipDistance(mfd)
{
  flippedAtoms = (bool*)malloc(numMolecules * sizeof(bool));
}

ClusterFlippedPoints::~ClusterFlippedPoints() {
  free(flippedAtoms);
}

void ClusterFlippedPoints::flip(int* molecules, vector<int> myFlipPoints) {
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  int iMolecule = molecules[myFlipPoints[1]];
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double* com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
  for (int j=1; j<myFlipPoints.size(); j++)
  {
    int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;
    int jMolecule = molecules[myFlipPoints[j]];
    box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom);

    for (int jAtom=jFirstAtom; jAtom<=jLastAtom; jAtom++) {
      double* r = box.getAtomPosition(jAtom);
      for (int i=0; i<3; i++) {
        r[i] = 2*com[i] - r[i];
      }
    }

  }
}

const double* ClusterFlippedPoints::getValues() {
  if (useCache && !cacheDirty) {
    return values;
  }
  cacheDirty = false;

  bool flipit = false;

  double minFlipDistance = 20;
  double minR2 = minFlipDistance*minFlipDistance;
  vector <vector <int>> myFlipPoints;

  for (int i=0; i<flipPoints.size(); i++)
  {
    flippedAtoms[i] = false;
    int idx0 = flipPoints[i][0];
    int idx1 = flipPoints[i][1];
    if (idx0 > idx1)
    {
      int t = idx0;
      idx0 = idx1;
      idx1 = t;
    }
    int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;

   box.getMoleculeInfo(idx0, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
   double* com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
    double com1[3] = {com[0],com[1],com[2]};
    box.getMoleculeInfo(idx1, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
    com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);

   double dr[3] = {com[0]-com1[0], com[1]-com1[1], com[2]-com1[2]};
   if (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] > minR2) {
     myFlipPoints.push_back(flipPoints[i]);
   }



  }


  const double* vsum = wrappedCluster.getValues();
  if (myFlipPoints.empty())
  {
    values = vsum;
    return vsum;
  }
  // for (int i=0; i<nValues; i++) {
  //   values[i] = foo[i];
  // }
  // if (!flipit) {
  //   return values;
  // }
  while (true) {
    bool didFlipTrue = false;
    for (int i=0; !didFlipTrue && i<myFlipPoints.size(); i++) {
      flippedAtoms[i] = !flippedAtoms[i];
      didFlipTrue = flippedAtoms[i];
      flip(atomList, myFlipPoints[i]);
    }
    if (!didFlipTrue) {
      // if we flipped every molecule from true to false, we must be done
      break;
    }
    const double* foo = wrappedCluster.getValues();
    vsum += foo;
  }

  values = vsum / pow(2, myFlipPoints.size());

  return values;
}
