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

void Species::init(AtomInfo& ai) {
  atomInfo = &ai;
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

// This could (almost) be static
double* Species::getMoleculeCOM(Box& box, int iFirstAtom, int iLastAtom) {
  com[0] = com[1] = com[2] = 0;
  double totMass = 0;
  double *r0 = box.getAtomPosition(iFirstAtom);
  for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    int iType = box.getAtomType(iAtom);
    double mass = atomInfo->getMass(iType);
    totMass += mass;
    if (iAtom==iFirstAtom) {
      for (int k=0; k<3; k++) com[k] = mass*ri[k];
    }
    else {
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = ri[k] - r0[k];
      box.nearestImage(dr);
      for (int k=0; k<3; k++) com[k] += mass*(r0[k] + dr[k]);
    }
  }
  for (int k=0; k<3; k++) com[k] /= totMass;
  return com;
}
