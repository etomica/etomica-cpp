#include <stdlib.h>
#include "species.h"

Species::Species(int na, int nat) : numAtoms(na), numAtomTypes(nat) {
  atomTypes = (int*)malloc(na*sizeof(int));
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
