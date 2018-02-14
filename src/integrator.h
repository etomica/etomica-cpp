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

class Integrator : public PotentialCallback {
  protected:
    PotentialMaster& potentialMaster;
    double temperature;
    double energy;
    int stepCount;
    vector<IntegratorListener*> listeners;
    vector<PotentialCallback*> selfPotentialCallbackVec;
  public:
    Integrator(PotentialMaster& potentialMaster);
    virtual ~Integrator() {}
    virtual void doStep() = 0;
    void doSteps(int steps);
    void setTemperature(double temperature);
    double getTemperature();
    void reset();
    double getPotentialEnergy();
    void addListener(IntegratorListener* listener);
    void allComputeFinished(double uTot, double virialTot, double** f);
};

class IntegratorMC : public Integrator {
  private:
    sfmt_t& sfmt;
    vector<MCMove*> moves;
    vector<double> moveProbabilities;
    double pMoveSum;
    MCMove* lastMove;
  public:
    IntegratorMC(PotentialMaster& potentialMaster, Random& random);
    IntegratorMC(PotentialMaster& potentialMaster, sfmt_t& random);
    ~IntegratorMC();
    void addMove(MCMove* move, double probability);
    virtual void doStep();
    void setTuning(bool doTuning);
    MCMove* getLastMove();
};

