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
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
};

class IntegratorMC : public Integrator {
  private:
    Random& random;
    vector<MCMove*> moves;
    vector<double> moveProbabilities;
    double pMoveSum;
    MCMove* lastMove;
  public:
    IntegratorMC(PotentialMaster& potentialMaster, Random& random);
    ~IntegratorMC();
    void addMove(MCMove* move, double probability);
    virtual void doStep();
    void setTuning(bool doTuning);
    MCMove* getLastMove();
};

#define THERMOSTAT_NONE 0
#define THERMOSTAT_NOSE_HOOVER 1
#define THERMOSTAT_ANDERSEN 2
#define THERMOSTAT_LANGEVIN 3

struct PotentialCallbackInfo {
  PotentialCallback* pcb;
  int interval;
  int countdown;
};

class IntegratorMD : public Integrator {
  protected:
    Random& random;
    Box& box;
    double** forces;
    double tStep;
    int thermostat;
    int nbrCheckInterval, nbrCheckCountdown;
    void randomizeVelocities(bool zeroMomentum);
    vector<struct PotentialCallbackInfo> allPotentialCallbacks;
  public:
    IntegratorMD(PotentialMaster& potentialMaster, Random& random, Box& box);
    ~IntegratorMD();
    void setTimeStep(double tStep);
    void setNbrCheckInterval(int interval);
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual void doStep();
    virtual void reset();
    void addPotentialCallback(PotentialCallback* callback, int interval=1);
};

