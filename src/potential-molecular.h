/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "potential.h"
#include "box.h"

class PotentialMolecular {
  public:
    PotentialMolecular() {}
    virtual ~PotentialMolecular() {}
    virtual double u(Box& box, int iFirstAtom, int jFirstAtom) = 0;
};

class PotentialMolecularAtomic : public PotentialMolecular {
  private:
    Potential& potential;
    SpeciesList& speciesList;

  public:
    PotentialMolecularAtomic(SpeciesList& sl, Potential& p);
    virtual ~PotentialMolecularAtomic() {}
    double u(Box& box, int iFirstAtom, int jFirstAtom);
};
