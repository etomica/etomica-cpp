/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "integrator.h"
#include "potential-master.h"
#include "potential-molecular.h"

class Cluster : public IntegratorListener {
  protected:
    const int nValues;
    const int numMolecules;
    bool useCache, cacheDirty, inTrial;
    double* values;
    double* oldValues;

  public:
    Cluster(int numMolecules, int numValues, bool cached);
    virtual ~Cluster();

    int numValues();
    void setCachingEnabled(bool enabled);
    virtual const double* getValues() = 0;
    void moveRejected(MCMove& move, double chi);
    void moveAccepted(MCMove& move, double chi);
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
    virtual double computeExtraFQ(int i);
  public:
    ClusterVirial(PotentialMasterVirial& potentialMaster, double temperature, int nDer, bool cached);
    virtual ~ClusterVirial();

    virtual const double* getValues();

};

class ClusterMBX : public ClusterVirial {
protected:
  PotentialMolecularMBnrg3body& p3body;
  virtual double computeExtraFQ(int i);
public:
  ClusterMBX(PotentialMasterVirial& potentialMaster, double temperature, int nDer, bool cached, PotentialMolecularMBnrg3body& p3body);
  virtual ~ClusterMBX();

};


class ClusterChain : public Cluster {
  protected:
    PotentialMasterVirial &potentialMaster;
    const double beta;
    int moleculePair[2];
    double chainFac, ringFac;

  public:
    ClusterChain(PotentialMasterVirial& potentialMaster, double temperature, double chainFac, double ringFac, bool cached);
    virtual ~ClusterChain();

    const double* getValues();
};

class ClusterFlipped : public Cluster {
  protected:
    Cluster& wrappedCluster;
    SpeciesList& speciesList;
    Box& box;
    bool* flippedAtoms;

    void flip(int iMolecule);

  public:
    ClusterFlipped(Cluster& cluster, SpeciesList& speciesList, Box& box, bool cached);
    virtual ~ClusterFlipped();

    const double* getValues();
};


class ClusterFlippedPoints : public Cluster {
protected:
  Cluster& wrappedCluster;
  SpeciesList& speciesList;
  Box& box;
  bool* flippedAtoms;
  vector<vector <int>> flipPoints;
  double minFlipDistance;

  void flip(vector<int> myFlipPoints);

public:
  ClusterFlippedPoints(Cluster& cluster, SpeciesList& speciesList, Box& box, bool cached, vector <vector<int>> fl, int mfd);
  virtual ~ClusterFlippedPoints();

  const double* getValues();
};
