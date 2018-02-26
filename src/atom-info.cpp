#include <stdlib.h>
#include <cstddef>

#include "atom-info.h"

AtomInfo::AtomInfo() : numAtomTypes(0), mass(nullptr) { }

AtomInfo::~AtomInfo() {
  if (mass) free(mass);
}

int AtomInfo::addAtomType(double m) {
  mass = (double*)realloc(mass, (numAtomTypes+1)*sizeof(double));
  mass[numAtomTypes] = m;
  numAtomTypes++;
  return numAtomTypes-1;
}

int AtomInfo::getNumTypes() const {
  return numAtomTypes;
}

double AtomInfo::getMass(int iType) const {
  return mass[iType];
}
