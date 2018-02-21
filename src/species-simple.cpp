#include "species.h"

SpeciesSimple::SpeciesSimple(int na, double m) : Species(na,1), mass(m) {
}

void SpeciesSimple::init(AtomInfo& atomInfo) {
  int myType = atomInfo.addAtomType(mass);
  for (int i=0; i<numAtoms; i++) atomTypes[i] = myType;
}
