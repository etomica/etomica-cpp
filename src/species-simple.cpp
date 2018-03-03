#include <algorithm>
#include "species.h"

SpeciesSimple::SpeciesSimple(int na, double m) : Species(na,1), mass(m) {
}

void SpeciesSimple::init(AtomInfo& atomInfo) {
  int myType = atomInfo.addAtomType(mass);
  for (int i=0; i<numAtoms; i++) {
    atomTypes[i] = myType;
  }
}

void SpeciesSimple::setAtomPosition(int iAtom, double* iPosition) {
  for (int k=0; k<3; k++) positions[iAtom][k] = iPosition[k];
}
