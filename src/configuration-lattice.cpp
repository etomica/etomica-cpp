/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <math.h>

#include "action.h"
#include "box.h"
#include "alloc2d.h"

void ConfigurationLattice::go(Box& box, int numBasisAtoms, double** basis, double* cellShape) {
  const double* boxSize = box.getBoxSize();
  int numMolecules = box.getTotalNumMolecules();
  int nCellsLeft = (numMolecules+(numBasisAtoms-1))/numBasisAtoms;
  int numCells[3] = {0,0,0};
  for (int dimLeft=3; dimLeft>0; dimLeft--) {
    double smin = 1e100;
    int dmin = 0;
    double product = 1.0;
    for (int idim = 0; idim < 3; idim++) {
      if (numCells[idim] > 0)
        continue;
      if (boxSize[idim]/cellShape[idim] < smin) {
        smin = boxSize[idim]/cellShape[idim];
        dmin = idim;
      }
      product *= boxSize[idim]/cellShape[idim];
    }
    // round off except for last dimension (then round up)
    if (dimLeft > 1) {
      numCells[dmin] = (int) round(boxSize[dmin]/cellShape[dmin] * pow((nCellsLeft / product), 1.0 / dimLeft));
    }
    else {
      numCells[dmin] = nCellsLeft;
    }
    if (numCells[dmin] == 0) {
      numCells[dmin] = 1;
    }
    nCellsLeft = (nCellsLeft + numCells[dmin] - 1) / numCells[dmin];
  }

  double myCellSize[3];
  for (int i=0; i<3; i++) myCellSize[i] = boxSize[i] / numCells[i];
  int ixyz[3];
  int iMolecule = 0;
  for (ixyz[0]=0; ixyz[0]<numCells[0]; ixyz[0]++) {
    for (ixyz[1]=0; ixyz[1]<numCells[1]; ixyz[1]++) {
      for (ixyz[2]=0; ixyz[2]<numCells[2]; ixyz[2]++) {
        for (int i=0; i<numBasisAtoms; i++) {
          if (iMolecule == numMolecules) break;
          double ri[3];
          for (int j=0; j<3; j++) {
            ri[j] = (basis[i][j] + ixyz[j] - 0.5*numCells[j]) * myCellSize[j];
          }
          int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
          box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
          Species* s = box.getSpeciesList().get(iSpecies);
          for (int jAtom=firstAtom; jAtom<=lastAtom; jAtom++) {
            double* jPos = s->getAtomPosition(jAtom-firstAtom);
            double* rj = box.getAtomPosition(jAtom);
            for (int k=0; k<3; k++) rj[k] = ri[k] + jPos[k];
            box.centralImage(rj);
          }
          iMolecule++;
        }
      }
    }
  }
}
