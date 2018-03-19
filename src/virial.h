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
    double newAlpha, newAlphaErr, alphaCor, alphaSpan;
    bool allDone, verbose, disposed;
    double alphaStats[4];
  public:
    VirialAlpha(IntegratorMC &refIntegrator, IntegratorMC &targetIntegrator, Cluster &refClusterRef, Cluster &refClusterTarget, Cluster &targetClusterRef, Cluster &targetClusterTarget);
    ~VirialAlpha();
    void setVerbose(bool newVerbose);
    void getNewAlpha(double &newAlpha, double &newAlphaErr, double &alphaCor);
    double* getAlphaStatistics();
    void setAlpha(double alphaCenter, double alphaSpan);
    void analyze(double &jBest);
    void runSteps(int steps);
    void run();
    bool getAllDone();
    void dispose();
    Average& getTargetAverage() {return targetAverage;}
    Average& getRefAverage() {return refAverage;}
};

class VirialProduction {
  protected:
    IntegratorMC &refIntegrator, &targetIntegrator;
    MeterVirialOverlap refMeter, targetMeter;
    AverageRatio refAverage, targetAverage;
    DataPump refPump, targetPump;
    double idealTargetFraction;
    double **refStats, **refBCStats, **refRatioStats;
    double **targetStats, **targetBCStats, **targetRatioStats;
    double alphaStats[2];
    double **fullStats;
    double refIntegral;
    bool disposed;
    long refSteps, targetSteps;
  public:
    VirialProduction(IntegratorMC &refIntegrator, IntegratorMC &targetIntegrator, Cluster &refClusterRef, Cluster &refClusterTarget, Cluster &targetClusterRef, Cluster &targetClusterTarget, double alpha, double refIntegral);
    ~VirialProduction();
    void dispose();
    void analyze();
    void printResults(const char **targetNames);
    void getResults();
    void runSteps(long numSteps);
    double** getFullStats();
    double* getAlphaStats();
    double** getRefStats();
    double** getTargetStats();
    double** getRefBCStats();
    double** getTargetBCStats();
    double** getRefRatioStats();
    double** getTargetRatioStats();
    AverageRatio& getTargetAverage() {return targetAverage;}
    AverageRatio& getRefAverage() {return refAverage;}
};
