#include <stdlib.h>
#include <math.h>
#include "species.h"
#include "alloc2d.h"

Species::Species(int na, int nat) : numAtoms(na), numAtomTypes(nat), atomTypes(nullptr), positions(nullptr) {
  if (na>0 && nat>0) setup(na, nat);
  axisAtoms[0][0] = axisAtoms[0][1] = axisAtoms[1][0] = axisAtoms[1][1] = -1;
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

void Species::getMoleculeOrientation(Box& box, int iFirstAtom, double* direction1, double* direction2) {
  if (axisAtoms[0][0] == -1) {
    direction1[0] = direction1[1] = direction1[2] = direction2[0] = direction2[1] = direction2[2] = 0;
    return;
  }
  double* r0 = box.getAtomPosition(iFirstAtom+axisAtoms[0][0]);
  double* r1 = box.getAtomPosition(iFirstAtom+axisAtoms[0][1]);
  for (int k=0; k<3; k++) direction1[k] = r1[k] - r0[k];
  box.nearestImage(direction1);
  double r2 = 0;
  for (int k=0; k<3; k++) r2 += direction1[k]*direction1[k];
  double r = sqrt(r2);
  for (int k=0; k<3; k++) direction1[k] /= r;
  if (axisAtoms[1][0] == -1) {
    direction2[0] = direction2[1] = direction2[2] = 0;
    return;
  }
  r0 = box.getAtomPosition(iFirstAtom+axisAtoms[1][0]);
  r1 = box.getAtomPosition(iFirstAtom+axisAtoms[1][1]);
  for (int k=0; k<3; k++) direction2[k] = r1[k] - r0[k];
  box.nearestImage(direction2);
  double dot = 0;
  for (int k=0; k<3; k++) dot += direction1[k]*direction2[k];
  for (int k=0; k<3; k++) direction2[k] -= dot*direction1[k];
  r2 = 0;
  for (int k=0; k<3; k++) r2 += direction2[k]*direction2[k];
  r = sqrt(r2);
  for (int k=0; k<3; k++) direction2[k] /= r;
}
