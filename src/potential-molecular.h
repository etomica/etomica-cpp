/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <json/json.hpp>
#include <potential/2b/x2b_A1B2_A1B2_v1x.h>
#include <potential/3b/mbnrg_3b_A1B2_A1B2_A1B2_deg4_v1.h>
#include <potential/dispersion/dispersion.h>
#include <potential/electrostatics/electrostatics.h>

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
  vector<double> grad;
  disp::Dispersion dispersion;
  elec::Electrostatics electrostatics;
  vector<double> chg;
  vector<double> chggrad;
  vector<double> pol;
  vector<double> polfac;
  double maxMin = 0.0;
  double minMax = 100;

public:
  PotentialMolecularMBnrg();
  virtual ~PotentialMolecularMBnrg();
  double u(Box& box, int iFirstAtom, int jFirstAtom);

};

class PotentialMolecularMBnrg3body : public PotentialMolecular
{
protected:
  double* xyz1{};
  double* xyz2{};
  double* xyz3{};
  mbnrg_A1B2_A1B2_A1B2_deg4::mbnrg_A1B2_A1B2_A1B2_deg4_v1 pot;
public:
  PotentialMolecularMBnrg3body();
  virtual ~PotentialMolecularMBnrg3body();
  double u(Box& box, int iFirstAtom, int jFirstAtom, int kFirstAtom);


};