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
