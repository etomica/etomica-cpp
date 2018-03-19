#pragma once

#include "move.h"
#include "box.h"
#include "random.h"
#include "SFMT.h"

class IntegratorListener {
  public:
    bool callFinished;
    bool callAccept, callReject;
    IntegratorListener();
    virtual ~IntegratorListener() {}
    virtual void stepFinished() {}
    virtual void moveAccepted(MCMove& move, double chi) {}
    virtual void moveRejected(MCMove& move, double chi) {}
};

class Integrator : public PotentialCallback {
  protected:
    PotentialMaster& potentialMaster;
    double temperature;
    double energy;
    long stepCount;
    vector<IntegratorListener*> listenersStepFinished;
    vector<PotentialCallback*> selfPotentialCallbackVec;
  public:
    Integrator(PotentialMaster& potentialMaster);
    virtual ~Integrator() {}
    virtual void doStep() = 0;
    void doSteps(int steps);
    long getStepCount();
    void setTemperature(double temperature);
    double getTemperature();
    void reset();
    double getPotentialEnergy();
    virtual void addListener(IntegratorListener* listener);
    virtual void removeListener(IntegratorListener* listener);
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
};

class IntegratorMC : public Integrator {
  private:
    Random& random;
    vector<MCMove*> moves;
    vector<double> moveProbabilities;
    double pMoveSum;
    MCMove* lastMove;
    vector<IntegratorListener*> listenersMoveAccepted, listenersMoveRejected;
  public:
    IntegratorMC(PotentialMaster& potentialMaster, Random& random);
    ~IntegratorMC();
    void addMove(MCMove* move, double probability);
    void removeMove(MCMove* move);
    virtual void doStep();
    void setTuning(bool doTuning);
    MCMove* getLastMove();
    virtual void addListener(IntegratorListener* listener);
    virtual void removeListener(IntegratorListener* listener);
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

