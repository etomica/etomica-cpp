#include <stdlib.h>
#include "species.h"
#include "alloc2d.h"

Species::Species(int na, int nat) : numAtoms(na), numAtomTypes(nat) {
  atomTypes = (int*)malloc(na*sizeof(int));
  positions = (double**)malloc2D(na, 3, sizeof(double));
  for (int i=0; i<na; i++) {
    for (int j=0; j<3; j++) positions[i][j] = 0;
  }
}

Species::~Species() {
  free(atomTypes);
}

int Species::getNumAtoms() {
  return numAtoms;
}

int* Species::getAtomTypes() {
  return atomTypes;
}

double* Species::getAtomPosition(int iAtom) {
  return positions[iAtom];
}
