#pragma once

#include "cluster.h"
#include "meter-virial.h"
#include "integrator.h"

class VirialAlpha {
  protected:
    long stepCount, nextCheck;
    IntegratorMC &refIntegrator, &targetIntegrator;
    MeterVirialOverlap &refMeter, &targetMeter;
    Average &refAverage, &targetAverage;
    double newAlpha, newAlphaErr, alphaCor;
    bool allDone, verbose;
  public:
    VirialAlpha(IntegratorMC &refIntegrator, IntegratorMC &targetIntegrator, MeterVirialOverlap &refMeter, MeterVirialOverlap &targetMeter, Average &refAverage, Average &targetAverage);
    ~VirialAlpha();
    void setVerbose(bool newVerbose);
    void getNewAlpha(double &newAlpha, double &newAlphaErr, double &alphaCor);
    void setAlpha(double alphaCenter, double alphaSpan);
    void analyze(double &jBest);
    void runSteps(int steps);
    void run();
};
