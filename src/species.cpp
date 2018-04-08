#include <stdlib.h>
#include "species.h"
#include "alloc2d.h"

Species::Species(int na, int nat) : numAtoms(na), numAtomTypes(nat), atomTypes(nullptr), positions(nullptr) {
  if (na>0 && nat>0) setup(na, nat);
}

void Species::setup(int na, int nat) {
  numAtoms = na;
  numAtomTypes = nat;
  atomTypes = (int*)realloc((void**)atomTypes, na*sizeof(int));
  positions = (double**)realloc2D((void**)positions, na, 3, sizeof(double));
  for (int i=0; i<na; i++) {
    for (int j=0; j<3; j++) positions[i][j] = 0;
  }
}

Species::~Species() {
  free(atomTypes);
  free2D((void**)positions);
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
