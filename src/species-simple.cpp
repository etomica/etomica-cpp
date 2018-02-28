#include <algorithm>
#include "species.h"

SpeciesSimple::SpeciesSimple(int na, double m) : Species(na,1), mass(m) {
}

void SpeciesSimple::init(AtomInfo& atomInfo) {
  int myType = atomInfo.addAtomType(mass);
  for (int i=0; i<numAtoms; i++) {
    atomTypes[i] = myType;
    for (int j=0; j<3; j++) positions[i][j] = 0;
  }
}

void SpeciesSimple::setAtomPosition(int iAtom, double* iPosition) {
  std::copy(positions[iAtom], positions[iAtom]+3, iPosition);
}
