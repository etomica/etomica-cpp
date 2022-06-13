/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <vector>
#include "potential-callback.h"
#include "random.h"
#include "SFMT.h"

using namespace std;

class AtomInfo;
class MCMove;
class PotentialMaster;
class Box;
class NeighborUpdateListener;

class IntegratorListener {
  public:
    bool callStepStarted, callStepFinished;
    bool callPreForce, callPostForce;
    bool callAccept, callReject;
    IntegratorListener();
    virtual ~IntegratorListener() {}
    virtual void stepStarted() {}
    virtual void stepFinished() {}
    virtual void preForce() {}
    virtual void postForce() {}
    virtual void moveAccepted(MCMove& move, double chi) {}
    virtual void moveRejected(MCMove& move, double chi) {}
};

class Integrator : public PotentialCallback {
  protected:
    PotentialMaster& potentialMaster;
    double temperature;
    double energy;
    long stepCount;
    vector<IntegratorListener*> listenersStepStarted;
    vector<IntegratorListener*> listenersStepFinished;
    vector<IntegratorListener*> listenersPreForce;
    vector<IntegratorListener*> listenersPostForce;
    vector<PotentialCallback*> selfPotentialCallbackVec;
  public:
    Integrator(PotentialMaster& potentialMaster);
    virtual ~Integrator() {}
    virtual void doStep() = 0;
    void doSteps(int steps);
    long getStepCount();
    void setTemperature(double temperature);
    double getTemperature();
    virtual void reset();
    double getPotentialEnergy();
    virtual void addListener(IntegratorListener* listener);
    virtual void removeListener(IntegratorListener* listener);
    virtual void allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor);
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
    AtomInfo& atomInfo;
    Random& random;
    Box& box;
    double** forces;
    double tStep;
    int thermostat;
    NeighborUpdateListener* neighborUpdateListener;
    vector<struct PotentialCallbackInfo> allPotentialCallbacks;
    double kineticEnergy;

    virtual void computeForces();
  public:
    IntegratorMD(AtomInfo& atomInfo, PotentialMaster& potentialMaster, Random& random, Box& box);
    virtual ~IntegratorMD();
    Box& getBox() {return box;}
    void setTimeStep(double tStep);
    double getTimeStep();
    double getTime() { return stepCount * tStep; }
    void setNbrCheckInterval(int interval);
    virtual void allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor);
    virtual void doStep() = 0;
    virtual void reset();
    void addPotentialCallback(PotentialCallback* callback, int interval=1);
    double getKineticEnergy();
    double** getForces();
    void randomizeVelocities(bool zeroMomentum);
    void randomizeVelocity(int iAtom);
};

class IntegratorMDVelocityVerlet : public IntegratorMD {
  protected:
  public:
    IntegratorMDVelocityVerlet(AtomInfo& atomInfo, PotentialMaster& potentialMaster, Random& random, Box& box);
    virtual ~IntegratorMDVelocityVerlet();
    virtual void doStep();
    virtual void reset();
};

class IntegratorNHC : public IntegratorMD {
  protected:
    const int numChains;
    double* q;
    double* eta;
    double* etaP;
  public:
    IntegratorNHC(AtomInfo& atomInfo, PotentialMaster& potentialMaster, Random& random, Box& box, int nChains, double tau);
    virtual ~IntegratorNHC();
    virtual void setTemperature(double newTermpature);
    virtual void doStep();
    virtual void propagatorU1(double dt);
    virtual void propagatorU2(double dt);
    virtual void propagatorU3(double dt);
    virtual void propagatorU4(double dt, int direction);
    virtual void reset();
};

#define ANDERSEN_DISABLED 0
#define ANDERSEN_SINGLE 1
#define ANDERSEN_FULL 2
#define ANDERSEN_PARTIAL 3
class AndersenThermostat : public IntegratorListener {
  private:
    IntegratorMD& integratorMD;
    Random& random;
    int mode;
    int interval;
    int intervalCountdown;
    double pRandomize;
    bool zeroMomentum;
  public:
    AndersenThermostat(IntegratorMD& integratorMD, Random& random);
    virtual ~AndersenThermostat();
    void setSingle(int interval);
    void setPartial(double pRandomize);
    void setFull(int interval, bool zeroMomentum);
    void stepFinished();
    double getTermperature();
    int getMode();
    void setTemperature(double temperature);
};

class IntegratorMulti : public Integrator {
  protected:
    vector<Integrator*> integrators;
    vector<long> steps;

  public:
    IntegratorMulti(PotentialMaster& potentialMasterDummy);
    virtual ~IntegratorMulti();
    void addIntegrator(Integrator* integrator, long steps);
    void setSteps(Integrator* integrator, long steps);
    Integrator* getIntegrator(int j);
    virtual void doStep() = 0;
    virtual void reset();
};
