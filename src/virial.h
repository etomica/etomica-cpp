#pragma once

#include "cluster.h"
#include "meter-virial.h"
#include "integrator.h"

class VirialAlpha {
  protected:
    long stepCount, nextCheck;
    IntegratorMC &refIntegrator, &targetIntegrator;
    MeterVirialOverlap refMeter, targetMeter;
    Average refAverage, targetAverage;
    DataPump refPump, targetPump;
    double newAlpha, newAlphaErr, alphaCor;
    bool allDone, verbose;
  public:
    VirialAlpha(IntegratorMC &refIntegrator, IntegratorMC &targetIntegrator, Cluster &refClusterRef, Cluster &refClusterTarget, Cluster &targetClusterRef, Cluster &targetClusterTarget);
    ~VirialAlpha();
    void setVerbose(bool newVerbose);
    void getNewAlpha(double &newAlpha, double &newAlphaErr, double &alphaCor);
    void setAlpha(double alphaCenter, double alphaSpan);
    void analyze(double &jBest);
    void runSteps(int steps);
    void run();
    Average& getTargetAverage() {return targetAverage;}
    Average& getRefAverage() {return refAverage;}
};

class VirialProduction {
  protected:
    long nextCheck;
    IntegratorMC &refIntegrator, &targetIntegrator;
    MeterVirialOverlap &refMeter, &targetMeter;
    AverageRatio &refAverage, &targetAverage;
    double idealTargetFraction;
    double **refStats, **refBCStats, **refRatioStats;
    double **targetStats, **targetBCStats, **targetRatioStats;
    double newAlpha, alphaErr;
    double *fullAvg, *fullErr;
    double refIntegral;
  public:
    VirialProduction(IntegratorMC &refIntegrator, IntegratorMC &targetIntegrator, MeterVirialOverlap &refMeter, MeterVirialOverlap &targetMeter, AverageRatio &refAverage, AverageRatio &targetAverage, double refIntegral);
    ~VirialProduction();
    void analyze();
    void printResults(const char **targetNames);
    void runSteps(long numSteps, int subSteps);
};
