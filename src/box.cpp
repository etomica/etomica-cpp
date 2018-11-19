/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "alloc2d.h"
#include "box.h"

Box::Box(SpeciesList &sl) : positions(nullptr), velocities(nullptr), knownNumSpecies(sl.size()), numAtomsBySpecies(nullptr), maxNumAtomsBySpecies(nullptr), speciesNumAtoms(nullptr), numMoleculesBySpecies(nullptr), maxNumMoleculesBySpecies(nullptr), firstAtom(nullptr), atomTypes(nullptr), speciesList(sl) {
  for (int i=0; i<3; i++) boxSize[i] = 0;

  int ss = knownNumSpecies;
  numAtomsBySpecies = new int[ss];
  maxNumAtomsBySpecies = new int[ss];
  speciesNumAtoms = new int[ss];
  numMoleculesBySpecies = new int[ss];
  maxNumMoleculesBySpecies = new int[ss];
  positions = (double***)malloc(ss*sizeof(double**));
  velocities = nullptr;
  firstAtom = (int**)malloc(ss*sizeof(int*));
  atomTypes = (int**)malloc(ss*sizeof(int*));
  for (int i=0; i<ss; i++) {
    positions[i] = nullptr;
    if (velocities) velocities[i] = nullptr;
    numMoleculesBySpecies[i] = numAtomsBySpecies[i] = maxNumMoleculesBySpecies[i] = 0;
    atomTypes[i] = firstAtom[i] = nullptr;
    speciesNumAtoms[i] = speciesList.get(i)->getNumAtoms();
  }
  periodic[0] = periodic[1] = periodic[2] = true;
}

Box::~Box() {
  if (knownNumSpecies==0) return;
  delete[] numAtomsBySpecies;
  delete[] numMoleculesBySpecies;
  delete[] maxNumMoleculesBySpecies;
  delete[] maxNumAtomsBySpecies;
  delete[] speciesNumAtoms;
  for (int i=0; i<knownNumSpecies; i++) {
    if (!positions[i]) continue;
    free2D((void**)positions[i]);
    if (velocities)free2D((void**)velocities[i]);
    free(firstAtom[i]);
    free(atomTypes[i]);
  }
  free(positions);
  free(velocities);
  free(firstAtom);
  free(atomTypes);
}

void Box::initCoordinates() {
  int dimLeft = 3;
  int numMolecules = getTotalNumMolecules();
  int nCellsLeft = (numMolecules+3)/4;
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
  double basisFCC[4][3];
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      basisFCC[i][j] = i==0?0:0.5;
    }
    if (i>0) basisFCC[i][i-1] = 0.0;
  }
  int ixyz[3];
  int iMolecule = 0;
  for (ixyz[0]=0; ixyz[0]<numCells[0]; ixyz[0]++) {
    for (ixyz[1]=0; ixyz[1]<numCells[1]; ixyz[1]++) {
      for (ixyz[2]=0; ixyz[2]<numCells[2]; ixyz[2]++) {
        for (int i=0; i<4; i++) {
          if (iMolecule == numMolecules) break;
          double ri[3];
          for (int j=0; j<3; j++) {
            ri[j] = (basisFCC[i][j] + ixyz[j] - 0.5*numCells[j]) * cellSize[j];
          }
          int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
          getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
          Species* s = speciesList.get(iSpecies);
          for (int jAtom=firstAtom; jAtom<=lastAtom; jAtom++) {
            double* jPos = s->getAtomPosition(jAtom-firstAtom);
            double* rj = getAtomPosition(jAtom);
            for (int k=0; k<3; k++) rj[k] = ri[k] + jPos[k];
            nearestImage(rj);
          }
          iMolecule++;
        }
      }
    }
  }
}

void Box::setNumMolecules(int iSpecies, int n) {
  if (iSpecies<0 || iSpecies>=knownNumSpecies) {
    fprintf(stderr, "unknown species id %d\n", iSpecies);
    abort();
  }
  Species* s = speciesList.get(iSpecies);
  int sna = s->getNumAtoms();
  int na = n*sna;
  if (n>maxNumMoleculesBySpecies[iSpecies]) {
    if (velocities) velocities[iSpecies] = (double**)realloc2D((void**)velocities[iSpecies], na, 3, sizeof(double));
    positions[iSpecies] = (double**)realloc2D((void**)positions[iSpecies], na, 3, sizeof(double));
    firstAtom[iSpecies] = (int*)realloc(firstAtom[iSpecies], n*sizeof(int));
    atomTypes[iSpecies] = (int*)realloc(atomTypes[iSpecies], na*sizeof(int));
    int* speciesAtomTypes = s->getAtomTypes();
    for (int i=numMoleculesBySpecies[iSpecies]; i<n; i++) {
      for (int j=i*sna; j<(i+1)*sna; j++) {
        atomTypes[iSpecies][j] = speciesAtomTypes[j-(i*sna)];
      }
    }
    maxNumMoleculesBySpecies[iSpecies] = n;
    maxNumAtomsBySpecies[iSpecies] = na;
  }
  int fa = 0;
  for (int jSpecies=0; jSpecies<knownNumSpecies; jSpecies++) {
    if (jSpecies<iSpecies) {
      fa += numAtomsBySpecies[jSpecies];
      continue;
    }
    int iStart = jSpecies==iSpecies ? numMoleculesBySpecies[iSpecies] : 0;
    fa += iStart;
    for (int i=iStart; i<n; i++) {
      firstAtom[iSpecies][i] = fa + i*sna;
    }
  }
  for (int iMolecule=numMoleculesBySpecies[iSpecies]; iMolecule<n; iMolecule++) {
    // place our new molecules at origin, using nominal conformation
    for (int j=0; j<sna; j++) {
      double* jPos = s->getAtomPosition(j);
      double* rj = positions[iSpecies][sna*iMolecule + j];
      rj[0] = jPos[0]; rj[1] = jPos[1]; rj[2] = jPos[2];
    }
  }
  numMoleculesBySpecies[iSpecies] = n;
  numAtomsBySpecies[iSpecies] = na;
}

double* Box::getAtomVelocity(int i) {
  int idx = i, iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && idx >= maxNumAtomsBySpecies[iSpecies]; iSpecies++) {
    idx -= maxNumAtomsBySpecies[iSpecies];
  }
#ifdef DEBUG
  if (idx>=maxNumAtomsBySpecies[iSpecies]) {
    printf("getAtomVecloities oops i %d is more atoms than I have\n", i);
    abort();
  }
#endif
  return velocities[0][i];
}

int Box::getGlobalMoleculeIndex(int iSpecies, int iMoleculeInSpecies) {
#ifdef DEBUG
  if (iSpecies>=knownNumSpecies) {
    printf("getGlobalMoleculeIndex oops i %d is more species than I have\n", iSpecies);
    abort();
  }
#endif
  int t = 0;
  for (int jSpecies=0; jSpecies<iSpecies; jSpecies++) {
    t += numMoleculesBySpecies[jSpecies];
  }
  return t + iMoleculeInSpecies;
}

void Box::getMoleculeInfoAtom(int iAtom, int &iMolecule, int &iSpecies, int &iFirstAtom) {
#ifdef DEBUG
  if (iAtom>=getNumAtoms()) {
    printf("getMoleculeInfoAtom oops i %d is more atoms than I have\n", iAtom);
    abort();
  }
#endif
  iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && iAtom >= numAtomsBySpecies[iSpecies]; iSpecies++) {
    iAtom -= numAtomsBySpecies[iSpecies];
  }
  iMolecule = iAtom/speciesNumAtoms[iSpecies];
  iFirstAtom = firstAtom[iSpecies][iMolecule];
}

void Box::getMoleculeInfo(int iMolecule, int &iSpecies, int &iMoleculeInSpecies, int &fa, int &la) {
  iMoleculeInSpecies = iMolecule;
  iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && iMoleculeInSpecies >= numMoleculesBySpecies[iSpecies]; iSpecies++) {
    iMoleculeInSpecies -= numMoleculesBySpecies[iSpecies];
  }
#ifdef DEBUG
  if (iMoleculeInSpecies>=numMoleculesBySpecies[iSpecies]) {
    printf("gFA oops i %d is more molecules than I have\n", iMolecule);
    abort();
  }
#endif
  fa = firstAtom[iSpecies][iMoleculeInSpecies];
  la = fa + speciesNumAtoms[iSpecies] - 1;
}

void Box::boxSizeUpdated() {
  for (int i=0; i<3; i++) boxHalf[i] = 0.5*boxSize[i];
}

void Box::nearestImage(double *dr) {
  for (int i=0; i<3; i++) {
    if (!periodic[i]) continue;
    while (dr[i] > boxHalf[i]) dr[i] -= boxSize[i];
    while (dr[i] < -boxHalf[i]) dr[i] += boxSize[i];
  }
}

const bool* Box::getPeriodic() {
  return periodic;
}

void Box::setPeriodic(const bool* newPeriodic) {
  for (int i=0; i<3; i++) periodic[i] = newPeriodic[i];
}

void Box::setBoxSize(double x, double y, double z) {
  boxSize[0] = x;
  boxSize[1] = y;
  boxSize[2] = z;
  boxSizeUpdated();
}

void Box::scaleBoxTo(double bx, double by, double bz) {
  double s[3] = {bx/boxSize[0]-1, by/boxSize[1]-1, bz/boxSize[2]-1};
  int nm = getTotalNumMolecules();
  // first unwrap and move molecules
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* center = species->getMoleculeCOM(*this, iFirstAtom, iLastAtom);
    for (int jAtom=iFirstAtom; jAtom<=iLastAtom; jAtom++) {
      double dr[3];
      double* rj = getAtomPosition(jAtom);
      for (int k=0; k<3; k++) dr[k] = rj[k]-center[k];
      nearestImage(dr);
      for (int k=0; k<3; k++) rj[k] = center[k]*s[k] + dr[k];
    }
  }
  // now change box size
  setBoxSize(bx, by, bz);
  int na = getNumAtoms();
  // now wrap atoms back inside box
  for (int iAtom=0; iAtom<na; iAtom++) {
    double* ri = getAtomPosition(iAtom);
    nearestImage(ri);
  }
}

void Box::enableVelocities() {
  if (velocities) {
    fprintf(stderr, "velocities alread enabled\n");
    return;
  }
  velocities = (double***)malloc(speciesList.size()*sizeof(double**));
  for (int i=0; i<speciesList.size(); i++) {
    int na = maxNumMoleculesBySpecies[i]*speciesList.get(i)->getNumAtoms();
    velocities[i] = (double**)malloc2D(na, 3, sizeof(double));
  }
}

