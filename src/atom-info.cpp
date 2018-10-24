/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
