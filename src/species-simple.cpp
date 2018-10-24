/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <algorithm>
#include "species.h"

SpeciesSimple::SpeciesSimple(int na, double m) : Species(na,1), mass(m) {
}

void SpeciesSimple::init(AtomInfo& ai) {
  Species::init(ai);
  int myType = ai.addAtomType(mass);
  for (int i=0; i<numAtoms; i++) {
    atomTypes[i] = myType;
  }
}

void SpeciesSimple::setAtomPosition(int iAtom, double* iPosition) {
  for (int k=0; k<3; k++) positions[iAtom][k] = iPosition[k];
}
