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
    const int nAtoms;
    Potential& potential;

  public:
    PotentialMolecularAtomic(int nAtoms, Potential& p);
    virtual ~PotentialMolecularAtomic() {}
    double u(Box& box, int iFirstAtom, int jFirstAtom);
};
