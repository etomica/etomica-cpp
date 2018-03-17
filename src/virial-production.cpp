#include <cmath>
#include "virial.h"
#include "alloc2d.h"

// perhaps just take Clusters and make Meter and Average internally

VirialProduction::VirialProduction(IntegratorMC &rIntegrator, IntegratorMC &tIntegrator, Cluster &refClusterRef, Cluster &refClusterTarget, Cluster &targetClusterRef, Cluster &targetClusterTarget, double alpha, double ri) : refIntegrator(rIntegrator), targetIntegrator(tIntegrator), refMeter(MeterVirialOverlap(refClusterRef, refClusterTarget, alpha, 0, 1)), targetMeter(MeterVirialOverlap(targetClusterTarget, targetClusterRef, 1/alpha, 0, 1)), refAverage(2, 1, 1000, true), targetAverage(targetClusterTarget.numValues()+1, 1, 1000, true), refPump(refMeter,1,&refAverage), targetPump(targetMeter,1,&targetAverage), idealTargetFraction(0.5), refIntegral(ri), refSteps(0), targetSteps(0) {
  int numTargets = targetAverage.getNumData();
  fullStats = (double**)malloc2D(numTargets, 2, sizeof(double));
  refIntegrator.addListener(&refPump);
  targetIntegrator.addListener(&targetPump);
}

VirialProduction::~VirialProduction() {
  dispose();
  free2D((void**)fullStats);
}

void VirialProduction::dispose() {
  if (disposed) return;
  refIntegrator.removeListener(&refPump);
  targetIntegrator.removeListener(&targetPump);
  disposed = true;
}

void VirialProduction::analyze() {
  refStats = refAverage.getStatistics();
  targetStats = targetAverage.getStatistics();
  if (std::isnan(targetStats[0][AVG_ERR])) {
    idealTargetFraction = 1;
    return;
  }
  if (std::isnan(refStats[0][AVG_ERR])) {
    idealTargetFraction = 0;
    return;
  }
  int numTargets = targetAverage.getNumData();
  double alpha = refMeter.getAlpha()[0];
  alphaStats[0] = refStats[1][AVG_AVG]/targetStats[numTargets-1][AVG_AVG]*alpha;
  alphaStats[1] = alpha*AverageRatio::ratioErr(refStats[1][AVG_AVG], refStats[1][AVG_ERR],
                               targetStats[numTargets-1][AVG_AVG], targetStats[numTargets-1][AVG_ERR], 0);

  refRatioStats = refAverage.getRatioStatistics();
  targetRatioStats = targetAverage.getRatioStatistics();
  refBCStats = refAverage.getBlockCorrelation();
  targetBCStats = targetAverage.getBlockCorrelation();
  for (int i=0; i<numTargets-1; i++) {
    fullStats[i][0] = refIntegral * targetRatioStats[i][AVG_AVG] / refRatioStats[0][AVG_AVG] * alpha;
    fullStats[i][1] = fabs(refIntegral) * AverageRatio::ratioErr(targetRatioStats[i][AVG_AVG], targetRatioStats[i][AVG_ERR], 
        refRatioStats[0][AVG_AVG], refRatioStats[0][AVG_ERR],0) * alpha;
    // also compute covariance between fullStats
    //printf("Bn(%d): %f %f\n", i, fullAvg[i], fullErr[i]);
    /*printf("ref: %f %f\n", refRatioStats[0][AVG_AVG], refRatioStats[0][AVG_ERR]);
    printf("  %f %f\n", refStats[0][AVG_AVG], refStats[0][AVG_ERR]);
    printf("  %f %f\n", refStats[1][AVG_AVG], refStats[1][AVG_ERR]);
    printf("tar: %f %f\n", targetRatioStats[i][AVG_AVG]*alpha, targetRatioStats[i][AVG_ERR]*alpha);
    printf("  %f %f\n", targetStats[0][AVG_AVG], targetStats[0][AVG_ERR]);
    printf("  %f %f\n", targetStats[1][AVG_AVG]/alpha, targetStats[1][AVG_ERR]/alpha);*/
  }

  double oldFrac = ((double)targetSteps)/(targetSteps + refSteps);
  double refErrorRatio = refRatioStats[0][AVG_ERR]/fabs(refRatioStats[0][AVG_AVG]);
  if (std::isnan(refErrorRatio) || refErrorRatio > 1) refErrorRatio = 1;
  double targetErrorRatio = targetRatioStats[0][AVG_ERR]/fabs(targetRatioStats[0][AVG_AVG]);
  if (std::isnan(targetErrorRatio) || targetErrorRatio > 1) targetErrorRatio = 1;
  idealTargetFraction = 1.0/(1 + refErrorRatio/targetErrorRatio * sqrt((1-oldFrac)/oldFrac));
}

void VirialProduction::printResults(const char **targetNames) {
  int numTargets = targetAverage.getNumData();
  printf("final reference step fraction: %5.4f\n", 1-idealTargetFraction);
  printf("actual reference step fraction: %5.4f\n", ((double)refSteps)/(refSteps+targetSteps));
  printf("alpha check:               % 22.15e  error: %12.5e\n", alphaStats[0], alphaStats[1]);
  printf("full average:              % 22.15e  error: %12.5e\n", fullStats[0][0], fullStats[0][1]);
  printf("reference ratio:           % 22.15e  error: %12.5e   cor: % 7.5f\n", refRatioStats[0][AVG_AVG], refRatioStats[0][AVG_ERR], refBCStats[0][1]);
  printf("reference average:         % 22.15e  error: %12.5e  acor: % 7.5f\n", refStats[0][AVG_AVG], refStats[0][AVG_ERR], refStats[0][AVG_ACOR]);
  printf("reference overlap average: % 22.15e  error: %12.5e  acor: % 7.5f\n", refStats[1][AVG_AVG], refStats[1][AVG_ERR], refStats[1][AVG_ACOR]);
  printf("target ratio:              % 22.15e  error: %12.5e   cor: % 7.5f\n", targetRatioStats[0][AVG_AVG], targetRatioStats[0][AVG_ERR], targetBCStats[0][1]);
  printf("target average:            % 22.15e  error: %12.5e  acor: % 7.5f\n", targetStats[0][AVG_AVG], targetStats[0][AVG_ERR], targetStats[0][AVG_ACOR]);
  printf("target overlap average:    % 22.15e  error: %12.5e  acor: % 7.5f\n", targetStats[numTargets-1][AVG_AVG], targetStats[numTargets-1][AVG_ERR], targetStats[numTargets-1][AVG_ACOR]);
  for (int i=1; i<numTargets-1; i++) {
    char name[28];
    if (targetNames && targetNames[i]) snprintf(name, 27, "%s average:", targetNames[i]);
    else snprintf(name, 27, "extra %d average:", i);
    printf("%-26s % 22.15e  error: %12.5e  acor: % 7.5f\n", name, targetStats[i][AVG_AVG], targetStats[i][AVG_ERR], targetStats[i][AVG_ACOR]);
    if (targetNames && targetNames[i]) snprintf(name, 27, "full %s average:", targetNames[i]);
    else snprintf(name, 27, "full extra %d average:", i);
    printf("%-26s % 22.15e  error: %12.5e\n", name, fullStats[i][0], fullStats[i][1]);
  }
}

double** VirialProduction::getFullStats() {return fullStats;}

double* VirialProduction::getAlphaStats() {return alphaStats;}

double** VirialProduction::getRefStats() { return refStats; }

double** VirialProduction::getTargetStats() {return targetStats;}

double** VirialProduction::getRefBCStats() {return refBCStats;}

double** VirialProduction::getTargetBCStats() {return targetBCStats;}

double** VirialProduction::getRefRatioStats() {return refRatioStats;}

double** VirialProduction::getTargetRatioStats() {return targetRatioStats;}

void VirialProduction::runSteps(long numSteps) {
  long totalSteps = refSteps + targetSteps;
  long subSteps = 100 + totalSteps/1000;
  if (subSteps > numSteps) subSteps = numSteps;
  long thisSteps = 0;
  while (thisSteps < numSteps) {
    bool runRef = true;
    if (totalSteps > 0) {
      double tFrac = ((double)targetSteps)/totalSteps;
      double idf = std::max(std::min(idealTargetFraction,0.99),0.01);
      runRef = tFrac > idf;
    }

    if (runRef) {
      refIntegrator.doSteps(subSteps);
      refSteps += subSteps;
    }
    else {
      targetIntegrator.doSteps(subSteps);
      targetSteps += subSteps;
    }

    analyze();
    totalSteps += subSteps;
    thisSteps += subSteps;

    subSteps = 100 + totalSteps/1000;
    if (subSteps > numSteps-thisSteps) subSteps = numSteps - thisSteps;
  }
}
