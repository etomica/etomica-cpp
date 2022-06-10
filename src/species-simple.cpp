/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <algorithm>
#include "species.h"

SpeciesSimple::SpeciesSimple(int na, int nat) : Species(na,nat) {
}

void SpeciesSimple::init(AtomInfo& ai) {
  Species::init(ai);
  int typeOffset = 0;
  for (int i=0; i<typeMass.size(); i++) {
    int myType = ai.addAtomType(typeMass[i]);
    typeOffset = myType-i;
  }
  int a = 0;
  for (int i=0; i<typeMass.size(); i++) {
    for (int j=0; j<numAtomsOfType[i]; j++) {
      atomTypes[a] = typeOffset + i;
      a++;
    }
  }
}

void SpeciesSimple::addAtomType(double mass, int atomsOfType) {
  typeMass.push_back(mass);
  numAtomsOfType.push_back(atomsOfType);
}

void SpeciesSimple::setAtomPosition(int iAtom, double x, double y, double z) {
  positions[iAtom][0] = x;
  positions[iAtom][1] = y;
  positions[iAtom][2] = z;
}
