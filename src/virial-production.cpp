#include "virial.h"

// perhaps just take Clusters and make Meter and Average internally

VirialProduction::VirialProduction(IntegratorMC &rIntegrator, IntegratorMC &tIntegrator, Cluster &refClusterRef, Cluster &refClusterTarget, Cluster &targetClusterRef, Cluster &targetClusterTarget, double alpha, double ri) : refIntegrator(rIntegrator), targetIntegrator(tIntegrator), refMeter(MeterVirialOverlap(refClusterRef, refClusterTarget, alpha, 0, 1)), targetMeter(MeterVirialOverlap(targetClusterTarget, targetClusterRef, alpha, 0, 1)), refAverage(2, 1, 1000, true), targetAverage(targetClusterTarget.numValues()+1, 1, 1000, true), refPump(refMeter,1,&refAverage), targetPump(targetMeter,1,&targetAverage), idealTargetFraction(0.5), refIntegral(ri) {
  int numTargets = targetAverage.getNumData();
  fullAvg = new double[numTargets];
  fullErr = new double[numTargets];
  refIntegrator.addListener(&refPump);
  targetIntegrator.addListener(&targetPump);
}

VirialProduction::~VirialProduction() {
  dispose();
  delete[] fullAvg;
  delete[] fullErr;
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
  newAlpha = refStats[1][AVG_AVG]/targetStats[numTargets-1][AVG_AVG]*alpha;
  alphaErr = alpha*AverageRatio::ratioErr(refStats[1][AVG_AVG], refStats[1][AVG_ERR],
                               targetStats[numTargets-1][AVG_AVG], targetStats[numTargets-1][AVG_ERR], 0);

  refBCStats = refAverage.getBlockCorrelation();
  targetBCStats = targetAverage.getBlockCovariance();
  refRatioStats = refAverage.getRatioStatistics();
  targetRatioStats = targetAverage.getRatioStatistics();
  for (int i=0; i<numTargets-1; i++) {
    fullAvg[i] = refIntegral * targetRatioStats[i][AVG_AVG] / refRatioStats[i][AVG_AVG] * alpha;
    fullErr[i] = fabs(refIntegral) * AverageRatio::ratioErr(targetRatioStats[i][AVG_AVG], targetRatioStats[i][AVG_ERR], 
        refRatioStats[0][AVG_AVG], refRatioStats[0][AVG_ERR],0) * alpha;
    // also compute covariance between fullAvg
    //printf("Bn(%d): %f %f\n", i, fullAvg[i], fullErr[i]);
    /*printf("ref: %f %f\n", refRatioStats[0][AVG_AVG], refRatioStats[0][AVG_ERR]);
    printf("  %f %f\n", refStats[0][AVG_AVG], refStats[0][AVG_ERR]);
    printf("  %f %f\n", refStats[1][AVG_AVG], refStats[1][AVG_ERR]);
    printf("tar: %f %f\n", targetRatioStats[i][AVG_AVG]*alpha, targetRatioStats[i][AVG_ERR]*alpha);
    printf("  %f %f\n", targetStats[0][AVG_AVG], targetStats[0][AVG_ERR]);
    printf("  %f %f\n", targetStats[1][AVG_AVG]/alpha, targetStats[1][AVG_ERR]/alpha);*/
  }

  double rSteps = refIntegrator.getStepCount();
  double tSteps = targetIntegrator.getStepCount();
  double oldFrac = tSteps/(tSteps + rSteps);
  double refErrorRatio = refRatioStats[0][AVG_ERR]/fabs(refRatioStats[0][AVG_AVG]);
  if (std::isnan(refErrorRatio) || refErrorRatio > 1) refErrorRatio = 1;
  double targetErrorRatio = targetRatioStats[0][AVG_ERR]/fabs(targetRatioStats[0][AVG_AVG]);
  if (std::isnan(targetErrorRatio) || targetErrorRatio > 1) targetErrorRatio = 1;
  idealTargetFraction = 1.0/(1 + refErrorRatio/targetErrorRatio * sqrt((1-oldFrac)/oldFrac));
}

void VirialProduction::printResults(const char **targetNames) {
  printf("final reference step fraction: %5.4f\n", 1-idealTargetFraction);
  double rSteps = refIntegrator.getStepCount();
  double tSteps = targetIntegrator.getStepCount();
  printf("actual reference step fraction: %5.4f\n", rSteps/(rSteps+tSteps));
  const char* name0 = targetNames && targetNames[0] ? targetNames[0] : "Full";
  printf("alpha check:               % 22.15e  error: %12.5e\n", newAlpha, alphaErr);
  printf("full average:              % 22.15e  error: %12.5e\n", fullAvg[0], fullErr[0]);
  printf("reference ratio:           % 22.15e  error: %12.5e   cor: % 7.5f\n", refRatioStats[0][AVG_AVG], refRatioStats[0][AVG_ERR], refBCStats[0][1]);
  printf("reference average:         % 22.15e  error: %12.5e  acor: % 7.5f\n", refStats[0][AVG_AVG], refStats[0][AVG_ERR], refStats[0][AVG_ACOR]);
  printf("reference overlap average: % 22.15e  error: %12.5e  acor: % 7.5f\n", refStats[1][AVG_AVG], refStats[1][AVG_ERR], refStats[1][AVG_ACOR]);
  printf("target ratio:              % 22.15e  error: %12.5e   cor: % 7.5f\n", targetRatioStats[0][AVG_AVG], targetRatioStats[0][AVG_ERR], targetBCStats[0][1]);
  printf("target average:            % 22.15e  error: %12.5e  acor: % 7.5f\n", targetStats[0][AVG_AVG], targetStats[0][AVG_ERR], targetStats[0][AVG_ACOR]);
  printf("target overlap average:    % 22.15e  error: %12.5e  acor: % 7.5f\n", targetStats[1][AVG_AVG], targetStats[1][AVG_ERR], targetStats[1][AVG_ACOR]);
  
}

void VirialProduction::runSteps(long numSteps, int subSteps) {
  for (int i=0; i<numSteps/subSteps; i++) {
    bool runRef = true;
    double rSteps = refIntegrator.getStepCount();
    double tSteps = targetIntegrator.getStepCount();
    if (rSteps + tSteps > 0) {
      double tFrac = tSteps/(rSteps + tSteps);
      runRef = tFrac > idealTargetFraction;
    }

    if (runRef) refIntegrator.doSteps(subSteps);
    else targetIntegrator.doSteps(subSteps);

    analyze();
  }
}
