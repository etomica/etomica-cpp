/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <potential/2b/x2b_A1B2_A1B2_v1x.h>

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

class PotentialMolecularMBnrg : public PotentialMolecular
{
protected:
  double* xyz1{};
  double* xyz2{};
  x2b_A1B2_A1B2_deg5::x2b_A1B2_A1B2_v1x pot;
public:
  PotentialMolecularMBnrg();
  virtual ~PotentialMolecularMBnrg();
  double u(Box& box, int iFirstAtom, int jFirstAtom);

};