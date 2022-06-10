/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "cluster.h"
#include "meter-virial.h"
#include "integrator.h"
#include "data-pump.h"
#include "data-sink.h"

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
    vector<double*> allAlphaStats;
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
    int getNumSavedStats();
    double* getSavedStats(int i);
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
    double **fullBCStats;
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
    double** getFullStats() {return fullStats;}
    double* getAlphaStats() {return alphaStats;}
    double** getRefStats() { return refStats; }
    double** getTargetStats() {return targetStats;}
    double** getRefBCStats() {return refBCStats;}
    double** getTargetBCStats() {return targetBCStats;}
    double** getRefRatioStats() {return refRatioStats;}
    double** getTargetRatioStats() {return targetRatioStats;}
    double** getFullBCStats() {return fullBCStats;}
    long getRefSteps() {return refSteps;}
    long getTargetSteps() {return targetSteps;}
    AverageRatio& getTargetAverage() {return targetAverage;}
    AverageRatio& getRefAverage() {return refAverage;}
};
