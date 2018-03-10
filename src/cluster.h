#pragma once

#include "potential-master.h"

class Cluster {
  protected:
    PotentialMasterVirial &potentialMaster;
    const int numMolecules;
    const double beta;
    const int nDer;
    double* values;
    double* oldValues;
    int** binomial;
    int moleculePair[2];

  public:
    Cluster(PotentialMasterVirial& potnetialMaster, double temperature, int nDer);
    ~Cluster();

    double* value();
    double* oldValue();
    void acceptNewValue();
    void rejectNewValue();
};
