#pragma once

#include "potential-master.h"

class Cluster {
  protected:
    PotentialMasterVirial &potentialMaster;
    const int numMolecules;
    const double beta;
    const int nDer;
    double* values;
    int** binomial;
    int moleculePair[2];

  public:
    Cluster(PotentialMasterVirial& potnetialMaster, double temperature, int nDer);
    ~Cluster() {}

    double* value();
};
