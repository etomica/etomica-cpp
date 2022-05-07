/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "move.h"
#include "cluster.h"
#include "rotation-matrix.h"

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

  public:

    MCMoveMoleculeDisplacementVirial(SpeciesList& speciesList, int iSpecies, Box& box, PotentialMaster& potentialMaster, Random& random, double stepSize, Cluster &cluster);
    ~MCMoveMoleculeDisplacementVirial();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
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

  public:
    MCMoveChainVirial(SpeciesList& speciesList, Box& box, PotentialMaster& potentialMaster, Random& random, double sigma);
    ~MCMoveChainVirial();

    virtual bool doTrial();
    virtual double getChi(double temperature) {return 1.0;}
    virtual void acceptNotify(){}
    virtual void rejectNotify() {fprintf(stderr, "no rejection for chain"); abort();}
    virtual double energyChange() {return 0;}
};

