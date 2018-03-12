#pragma once

#include "potential-master.h"

class Cluster {
  protected:
    const int nValues;
    const int numMolecules;
    bool useCache, cacheDirty;
    double* values;
    double* oldValues;

  public:
    Cluster(int numMolecules, int numValues, bool cached);
    virtual ~Cluster();

    int numValues();
    void setCachingEnabled(bool enabled);
    virtual const double* getValues() = 0;
    void trialNotify();
    void trialRejected();
};

class ClusterVirial : public Cluster {
  protected:
    PotentialMasterVirial &potentialMaster;
    const double beta;
    const int nDer;
    int** binomial;
    int moleculePair[2];
    double prefac;

  public:
    ClusterVirial(PotentialMasterVirial& potentialMaster, double temperature, int nDer, bool cached);
    virtual ~ClusterVirial();

    const double* getValues();
};

class ClusterChain : public Cluster {
  protected:
    PotentialMasterVirial &potentialMaster;
    const double beta;
    int moleculePair[2];
    double chainFac, ringFac;

  public:
    ClusterChain(PotentialMasterVirial& potentialMaster, double temperature, double chainFac, double ringFac);
    virtual ~ClusterChain();

    const double* getValues();
};
