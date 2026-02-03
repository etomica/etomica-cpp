/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <utility>

#include "move.h"
#include "cluster.h"
#include "rotation-matrix.h"
#include "vector.h"

class MCMoveDisplacementVirial : public MCMove {
  private:
    Cluster &cluster;
    double rOld[3];
    double wOld, wNew;
    int iAtom;

  public:

    MCMoveDisplacementVirial(Box& box, PotentialMaster& potentialMaster, Random& random, double stepSize, Cluster &cluster);
    ~MCMoveDisplacementVirial();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
};

class MCMoveMoleculeDisplacementVirial : public MCMove {
  private:
    Cluster &cluster;
    double **rOld;
    double wOld, wNew;
    int iMolecule;
    int iSpecies;
    double pisum[91], piHist[91], histogram[91];
    long hcount[91];
    void addToHistogram(double pi);

  public:

    MCMoveMoleculeDisplacementVirial(SpeciesList& speciesList, int iSpecies, Box& box, PotentialMaster& potentialMaster, Random& random, double stepSize, Cluster &cluster);
    ~MCMoveMoleculeDisplacementVirial();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
    double* getHistogram();
    double* getHistogramPi();
};

class MCMoveMoleculeRotateVirial : public MCMove {
  private:
    SpeciesList& speciesList;
    Cluster &cluster;
    RotationMatrix mat;
    double **rOld;
    double wOld, wNew;
    int iMolecule;
    int iSpecies;
    int nTimesPointedX = 0;
    int total = 0;
    bool first = true;

  public:

    MCMoveMoleculeRotateVirial(SpeciesList& speciesList, int iSpecies, Box& box, PotentialMaster& potentialMaster, Random& random, double stepSize, Cluster &cluster);
    ~MCMoveMoleculeRotateVirial();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
};

class MCMoveChainVirial : public MCMove {
  private:
    SpeciesList& speciesList;
    const double sigma;
    double histogram[61];
    long hcount[61];

  public:
    MCMoveChainVirial(SpeciesList& speciesList, Box& box, PotentialMaster& potentialMaster, Random& random, double sigma);
    ~MCMoveChainVirial();

    virtual bool doTrial();
    virtual double getChi(double temperature) {return 1.0;}
    virtual void acceptNotify(){}
    virtual void rejectNotify() {fprintf(stderr, "no rejection for chain"); abort();}
    virtual double energyChange() {return 0;}
    double* getHistogram();
};

class MCMoveClusterAngleGeneral : public MCMove
{
  private:
    SpeciesList& speciesList;
    bool oneSide;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    vector<int> modified;
    int b = 0;
    Cluster &cluster;
    RotationMatrix mat;
    int iAtomFirst{}, iAtomLast{};
    int numOldPositions{};
    double **oldPositions{};
    int counter =0;
    double sum = 0;

    vector<int *> triplets;
    vector<vector <int>> bonding;

    int iMolecule;
    int mySpecies = -1;
    int* constraintMap{};
    void transformBondedAtoms(int index, double* shift);
    void transform(int index, double* shift);


  public:
    int idx = -1;
    MCMoveClusterAngleGeneral(Box &b, PotentialMaster &p, bool oneSide, Random &r, vector<vector <int>> bnd, vector<int *> t, SpeciesList& sl,
                              double stepSize, Cluster &cluster);

    ~MCMoveClusterAngleGeneral();

    // void setBox(Box& p);
    void setConstraintMap(vector<int> newConstraintMap) {constraintMap = newConstraintMap.data();}
    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
};

class MCMoveClusterStretch : public MCMove
{
private:
  SpeciesList& speciesList;
  double uOld = 0;
  double uNew = 0;
  double wOld = 0;
  double wNew = 0;
  vector<int> modified;
  int b = 0;
  Cluster &cluster;
  int iAtomFirst{}, iAtomLast{};
  int numOldPositions{};
  double **oldPositions{};
  int counter =0;
  double sum = 0;
  double step = 0;
  vector<vector <int>> bonding;
  int iMolecule;
  int mySpecies = -1;
  double transformBondedAtoms(double* dr, int index, double* shift);
  double transform(double* dr, int index, double* shift);


public:
  int idx = -1;
  MCMoveClusterStretch(Box &b, PotentialMaster &p, Random &r, vector<vector <int>> bnd, SpeciesList& sl,
                            double stepSize, Cluster &cluster);

  ~MCMoveClusterStretch();

  // void setBox(Box& p);
  virtual bool doTrial();
  virtual double getChi(double temperature);
  virtual void acceptNotify();
  virtual void rejectNotify();
  virtual double energyChange();
};