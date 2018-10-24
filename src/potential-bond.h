/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "potential.h"

class PotentialBondHarmonic : public Potential {
  protected:
    double r0;
    double epsilon;

  public:
    PotentialBondHarmonic(double r0, double eps);
    virtual ~PotentialBondHarmonic();
    virtual double ur(double r);
    virtual double u(double r2);
    virtual double du(double r2);
    virtual double d2u(double r2);
    virtual void u012(double r2, double &u, double &du, double &d2u);
};
