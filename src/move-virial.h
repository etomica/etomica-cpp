#pragma once

#include "move.h"
#include "cluster.h"

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

class MCMoveChainVirial : public MCMove {
  private:
    const double sigma;

  public:
    MCMoveChainVirial(Box& box, PotentialMaster& potentialMaster, Random& random, double sigma);
    ~MCMoveChainVirial();

    virtual bool doTrial();
    virtual double getChi(double temperature) {return 1.0;}
    virtual void acceptNotify(){}
    virtual void rejectNotify() {fprintf(stderr, "no rejection for chain"); abort();}
    virtual double energyChange() {return 0;}
};

