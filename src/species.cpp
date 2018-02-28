#include <stdlib.h>
#include "species.h"
#include "alloc2d.h"

Species::Species(int na, int nat) : numAtoms(na), numAtomTypes(nat) {
  atomTypes = (int*)malloc(na*sizeof(int));
  positions = (double**)malloc2D(na, 3, sizeof(double));
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
