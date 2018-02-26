#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "alloc2d.h"
#include "box.h"

Box::Box(SpeciesList &sl) : positions(nullptr), velocities(nullptr), knownNumSpecies(sl.size()), numAtomsBySpecies(nullptr), numMoleculesBySpecies(nullptr), maxNumMoleculesBySpecies(nullptr), firstAtom(nullptr), moleculeIdx(nullptr), atomTypes(nullptr), speciesList(sl) {
  for (int i=0; i<3; i++) boxSize[i] = 0;

  int ss = knownNumSpecies;
  numAtomsBySpecies = (int*)realloc(numAtomsBySpecies, ss*sizeof(int));
  numMoleculesBySpecies = (int*)realloc(numMoleculesBySpecies, ss*sizeof(int));
  maxNumMoleculesBySpecies = (int*)realloc(maxNumMoleculesBySpecies, ss*sizeof(int));
  positions = (double***)realloc(positions, ss*sizeof(double**));
  if (velocities) velocities = (double***)realloc(velocities, ss*sizeof(double**));
  firstAtom = (int**)realloc(firstAtom, ss*sizeof(int*));
  moleculeIdx = (int**)realloc(moleculeIdx, ss*sizeof(int*));
  atomTypes = (int**)realloc(atomTypes, ss*sizeof(int*));
  if (velocities) velocities = (double***)realloc(velocities, ss*sizeof(double**));
  for (int i=0; i<ss; i++) {
    positions[i] = nullptr;
    if (velocities) velocities[i] = nullptr;
    numMoleculesBySpecies[i] = numAtomsBySpecies[i] = maxNumMoleculesBySpecies[i] = 0;
    atomTypes[i] = firstAtom[i] = moleculeIdx[i] = nullptr;
  }
}

Box::~Box() {
  if (knownNumSpecies==0) return;
  free(numAtomsBySpecies);
  free(numMoleculesBySpecies);
  free(maxNumMoleculesBySpecies);
  for (int i=0; i<knownNumSpecies; i++) {
    free2D((void**)positions[i]);
    if (velocities)free2D((void**)velocities[i]);
    free(firstAtom[i]);
    free(moleculeIdx[i]);
  }
  free(positions);
  free(velocities);
  free(firstAtom);
  free(moleculeIdx);
  free(atomTypes);
}

int Box::getNumMolecules(int iSpecies) {
  return numMoleculesBySpecies[iSpecies];
}

int Box::getTotalNumMolecules() {
  int s = 0;
  for (int i=0; i<knownNumSpecies; i++) s += numMoleculesBySpecies[i];
  return s;
}

int Box::getNumAtoms() {
  int s = 0;
  for (int i=0; i<knownNumSpecies; i++) s += numAtomsBySpecies[i];
  return s;
}

void Box::initCoordinates() {
  int dimLeft = 3;
  int numAtoms = getNumAtoms();
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
          if (iAtom == numAtoms) break;
          double* ri = getAtomPosition(iAtom);
          for (int j=0; j<3; j++) {
            ri[j] = (basisFCC[i][j] + ixyz[j] - 0.5*numCells[j]) * cellSize[j];
          }
          iAtom++;
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
    if (iSpecies==0) {
      int newTotalAtoms = getNumAtoms() + (n-numMoleculesBySpecies[iSpecies])*sna;
      positions[0] = (double**)realloc2D0((void**)positions[0], na, newTotalAtoms, 3, sizeof(double));
      if (velocities) velocities[0] = (double**)realloc2D0((void**)velocities[0], na, newTotalAtoms, 3, sizeof(double));
      // positions[0] will properly point to the raw positions for species 0.  it has memory allocated
      // to point to raw positions for other species, but those pointers have not been assigned.
      // do that now.
      int jAtom = na;
      for (int jSpecies=1; jSpecies<knownNumSpecies; jSpecies++) {
        for (int j=0; j<numAtomsBySpecies[jSpecies]; j++) {
          positions[0][jAtom] = positions[jSpecies][j];
          velocities[0][jAtom] = velocities[jSpecies][j];
          jAtom++;
        }
      }
    }
    else {
      if (velocities) velocities[iSpecies] = (double**)realloc2D((void**)velocities[iSpecies], na, 3, sizeof(double));
      positions[iSpecies] = (double**)realloc2D((void**)positions[iSpecies], na, 3, sizeof(double));
      // positions[0] needs to be updated to point to our new positions as well
      int jAtom = 0;
      for (int jSpecies=0; jSpecies<knownNumSpecies; jSpecies++) {
        if (jSpecies<iSpecies) {
          jAtom += numAtomsBySpecies[jSpecies];
          continue;
        }
        for (int j=0; j<numAtomsBySpecies[jSpecies]; j++) {
          positions[0][jAtom] = positions[jSpecies][j];
          velocities[0][jAtom] = velocities[jSpecies][j];
          jAtom++;
        }
      }
    }
    firstAtom[iSpecies] = (int*)realloc(firstAtom[iSpecies], n*sizeof(int));
    moleculeIdx[iSpecies] = (int*)realloc(moleculeIdx[iSpecies], n*sizeof(int));
    atomTypes[iSpecies] = (int*)realloc(atomTypes[iSpecies], n*sizeof(int));
    int* speciesAtomTypes = s->getAtomTypes();
    for (int i=numMoleculesBySpecies[iSpecies]; i<n; i++) {
      firstAtom[iSpecies][i] = i*sna;
      for (int j=i*sna; j<(i+1)*sna; j++) {
        moleculeIdx[iSpecies][j] = i;
        atomTypes[iSpecies][j] = speciesAtomTypes[j-(i*sna)];
      }
    }
    maxNumMoleculesBySpecies[iSpecies] = n;
  }
  numMoleculesBySpecies[iSpecies] = n;
  numAtomsBySpecies[iSpecies] = na;
}

double* Box::getAtomPosition(int i) {
#ifdef DEBUG
  int idx = i, iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && idx > numAtomsBySpecies[iSpecies]; iSpecies++) {
    idx -= numAtomsBySpecies[iSpecies];
  }
  if (idx>=numAtomsBySpecies[iSpecies]) {
    printf("gAP oops i %d is more atoms than I have\n", i);
    abort();
  }
#endif
  return positions[0][i];
}

double* Box::getAtomVelocity(int i) {
#ifdef DEBUG
  int idx = i, iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && idx > numAtomsBySpecies[iSpecies]; iSpecies++) {
    idx -= numAtomsBySpecies[iSpecies];
  }
  if (idx>=numAtomsBySpecies[iSpecies]) {
    printf("gAP oops i %d is more atoms than I have\n", i);
    abort();
  }
#endif
  return velocities[0][i];
}

int Box::getAtomType(int i) {
  int idx = i, iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && idx > numAtomsBySpecies[iSpecies]; iSpecies++) {
    idx -= numAtomsBySpecies[iSpecies];
  }
#ifdef DEBUG
  if (idx>=numAtomsBySpecies[iSpecies]) {
    printf("gAT oops i %d is more atoms than I have\n", i);
    abort();
  }
#endif
  return atomTypes[iSpecies][idx];
}

int Box::getMolecule(int iAtom) {
  int idx = iAtom, iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && idx > numAtomsBySpecies[iSpecies]; iSpecies++) {
    idx -= numAtomsBySpecies[iSpecies];
  }
#ifdef DEBUG
  if (idx>=numAtomsBySpecies[iSpecies]) {
    printf("gAT oops i %d is more atoms than I have\n", i);
    abort();
  }
#endif
  return firstAtom[iSpecies][idx];
}

int Box::getFirstAtom(int iSpecies, int iMoleculeInSpecies) {
  return firstAtom[iSpecies][iMoleculeInSpecies];
}

void Box::getMoleculeInfo(int iMolecule, int &iSpecies, int &fa, int &la) {
  int idx = iMolecule;
  iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && idx > numMoleculesBySpecies[iSpecies]; iSpecies++) {
    idx -= numMoleculesBySpecies[iSpecies];
  }
#ifdef DEBUG
  if (idx>=numMoleculesBySpecies[iSpecies]) {
    printf("gFA oops i %d is more molecules than I have\n", i);
    abort();
  }
#endif
  fa = firstAtom[iSpecies][idx];
  la = fa + speciesList.get(iSpecies)->getNumAtoms() - 1;
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

void Box::enableVelocities() {
  if (velocities) {
    fprintf(stderr, "velocities alread enabled\n");
    return;
  }
  velocities = (double***)realloc(velocities, speciesList.size()*sizeof(double**));
  for (int i=0; i<speciesList.size(); i++) {
    int na = maxNumMoleculesBySpecies[i]*speciesList.get(i)->getNumAtoms();
    velocities[i] = (double**)realloc2D((void**)velocities[i], na, 3, sizeof(double));
  }
}

