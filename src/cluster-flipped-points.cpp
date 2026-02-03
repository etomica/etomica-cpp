/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "cluster.h"

ClusterFlippedPoints::ClusterFlippedPoints(Cluster& c, SpeciesList& sl, Box& b, bool cached, vector <vector<int>> fl, int mfd) :
  Cluster(b.getTotalNumMolecules(), c.numValues(), cached), wrappedCluster(c), speciesList(sl), box(b),
   flipPoints(fl), minFlipDistance(mfd)
{
  flippedAtoms = (bool*)malloc(numMolecules * sizeof(bool));
}

ClusterFlippedPoints::~ClusterFlippedPoints() {
  free(flippedAtoms);
}

void ClusterFlippedPoints::flip(vector<int> myFlipPoints) {
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  int iMolecule = myFlipPoints[1];
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double* com = speciesList.get(iSpecies)->getMoleculeCOM(box, firstAtom, lastAtom);
  for (int j=1; j<myFlipPoints.size(); j++)
  {
    int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;
    int jMolecule = myFlipPoints[j];
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
    // printf("dr: %f %f %f \n", dr[0], dr[1], dr[2]);
   if (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] > minR2) {
     myFlipPoints.push_back(flipPoints[i]);
   }



  }


  const double* vsum = wrappedCluster.getValues();
  for (int i = 0; i < wrappedCluster.numValues(); i++)
  {
    values[i] = vsum[i];
  }

  if (myFlipPoints.empty())
  {
    return values;
  }
  // values[0] = 0;
  // return values;
  // printf("vsum: %f \n", vsum[0]);
  while (true) {
    bool didFlipTrue = false;
    for (int i=0; !didFlipTrue && i<myFlipPoints.size(); i++) {
      flippedAtoms[i] = !flippedAtoms[i];
      didFlipTrue = flippedAtoms[i];
      flip(myFlipPoints[i]);
    }
    if (!didFlipTrue) {
      // if we flipped every molecule from true to false, we must be done
      break;
    }
    vsum = wrappedCluster.getValues();
    // printf("%f \n", vsum[0]);
    for (int i = 0; i < wrappedCluster.numValues(); i++)
    {
      values[i] += vsum[i];
    }
  }
  for (int i = 0; i < wrappedCluster.numValues(); i++)
  {
    values[i] /= pow(2, myFlipPoints.size());

  }
  // printf("values: %f \n", values[0]);

  return values;
}
