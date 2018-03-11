#pragma once

#include "potential-master.h"

class Cluster {
  protected:
    PotentialMasterVirial &potentialMaster;
    const int numMolecules;
    const double beta;
    const int nDer;
    bool useCache, cacheDirty;
    double* values;
    double* oldValues;
    int** binomial;
    int moleculePair[2];

  public:
    Cluster(PotentialMasterVirial& potnetialMaster, double temperature, int nDer, bool cached);
    ~Cluster();

    int numValues();
    void setCachingEnabled(bool enabled);
    const double* getValues();
    void trialNotify();
    void trialRejected();
};
