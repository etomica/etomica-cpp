#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "alloc2d.h"
#include "box.h"

Box::Box() {
  for (int i=0; i<3; i++) boxSize[i] = 0;
  coords = nullptr;
  numAtoms = 0;
  maxNumAtoms = 0;
}

Box::~Box() {
  free2D((void**)coords);
}

int Box::getNumAtoms() {
  return numAtoms;
}

void Box::initCoordinates() {
  int dimLeft = 3;
  int nCellsLeft = (numAtoms+3)/4;
  int numCells[3] = {0,0,0};
  while (dimLeft > 0) {
    double smin = 1e100;
    int dmin = 0;
    double product = 1.0;
    for (int idim = 0; idim < 3; idim++) {
      if (numCells[idim] > 0)
        continue;
      if (boxSize[idim] < smin) {
        smin = boxSize[idim];
        dmin = idim;
      }
      product *= boxSize[idim];
    }
    // round off except for last dimension (then round up)
    if (dimLeft > 1) {
      numCells[dmin] = (int) round(boxSize[dmin] * pow((nCellsLeft / product), 1.0 / dimLeft));
    }
    else {
      numCells[dmin] = nCellsLeft;
    }
    if (numCells[dmin] == 0) {
      numCells[dmin] = 1;
    }
    nCellsLeft = (nCellsLeft + numCells[dmin] - 1) / numCells[dmin];
    dimLeft--;
  }

  double cellSize[3];
  for (int i=0; i<3; i++) cellSize[i] = boxSize[i] / numCells[i];
  int iAtom = 0;
  double basisFCC[4][3];
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      basisFCC[i][j] = i==0?0:0.5;
    }
    if (i>0) basisFCC[i][i-1] = 0.0;
  }
  int ixyz[3];
  for (ixyz[0]=0; ixyz[0]<numCells[0]; ixyz[0]++) {
    for (ixyz[1]=0; ixyz[1]<numCells[1]; ixyz[1]++) {
      for (ixyz[2]=0; ixyz[2]<numCells[2]; ixyz[2]++) {
        for (int i=0; i<4; i++) {
          for (int j=0; j<3; j++) {
            coords[iAtom][j] = (basisFCC[i][j] + ixyz[j] - 0.5*numCells[j]) * cellSize[j];
          }
          iAtom++;
          if (iAtom == numAtoms) return;
        }
      }
    }
  }
}

void Box::setNumAtoms(int n) {
  if (n>maxNumAtoms) {
    coords = (double**)realloc2D((void**)coords, n, 3, sizeof(double));
    maxNumAtoms = n;
  }
  numAtoms = n;
}

double* Box::getAtomPosition(int i) {
  /*if (numAtoms != 500) {
    printf("gAP oops %d %d (%d)\n", i, numAtoms, maxNumAtoms);
    abort();
  }
  if (i>=numAtoms) {
    printf("gAP oops %d > %d (%d)\n", i, numAtoms, maxNumAtoms);
    abort();
  }*/
  return coords[i];
}

void Box::boxSizeUpdated() {
  for (int i=0; i<3; i++) boxHalf[i] = 0.5*boxSize[i];
}

void Box::nearestImage(double *dr) {
  for (int i=0; i<3; i++) {
    while (dr[i] > boxHalf[i]) dr[i] -= boxSize[i];
    while (dr[i] < -boxHalf[i]) dr[i] += boxSize[i];
  }
}

void Box::setBoxSize(double x, double y, double z) {
  boxSize[0] = x;
  boxSize[1] = y;
  boxSize[2] = z;
  boxSizeUpdated();
}
