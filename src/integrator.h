#pragma once

#include "move.h"
#include "box.h"
#include "random.h"
#include "SFMT.h"

class IntegratorListener {
  public:
    IntegratorListener();
    virtual ~IntegratorListener() {}
    virtual void stepFinished() = 0;
};

class IntegratorMC : public PotentialCallback {
  private:
    PotentialMaster& potentialMaster;
    sfmt_t& sfmt;
    double temperature;
    double energy;
    int stepCount;
    vector<MCMove*> moves;
    vector<double> moveProbabilities;
    vector<IntegratorListener*> listeners;
    vector<PotentialCallback*> selfPotentialCallbackVec;
    double pMoveSum;
  public:
    IntegratorMC(PotentialMaster& potentialMaster, Random& random);
    IntegratorMC(PotentialMaster& potentialMaster, sfmt_t& random);
    ~IntegratorMC();
    void addMove(MCMove* move, double probability);
    void doStep();
    void doSteps(int steps);
    void setTemperature(double temperature);
    double getTemperature();
    void reset();
    double getPotentialEnergy();
    void addListener(IntegratorListener* listener);
    void allComputeFinished(double uTot, double virialTot, double** f);
    void setTuning(bool doTuning);
};

