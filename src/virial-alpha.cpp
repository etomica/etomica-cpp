#include "virial.h"

// perhaps just take Clusters and make Meter and Average internally

VirialAlpha::VirialAlpha(IntegratorMC &rIntegrator, IntegratorMC &tIntegrator, Cluster &refClusterRef, Cluster &refClusterTarget, Cluster &targetClusterRef, Cluster &targetClusterTarget) : stepCount(0), nextCheck(1000), refIntegrator(rIntegrator), targetIntegrator(tIntegrator), refMeter(MeterVirialOverlap(refClusterRef, refClusterTarget, 1, 5, 10)), targetMeter(MeterVirialOverlap(targetClusterTarget, targetClusterRef, 1, -5, 10)), refAverage(10, 1, 1000, false), targetAverage(10, 1, 100, false), refPump(refMeter,1,&refAverage), targetPump(targetMeter,1,&targetAverage), allDone(false), verbose(false) {
  refIntegrator.addListener(&refPump);
  targetIntegrator.addListener(&targetPump);
}

VirialAlpha::~VirialAlpha() {
  refIntegrator.removeListener(&refPump);
  targetIntegrator.removeListener(&targetPump);
}

void VirialAlpha::setVerbose(bool newVerbose) {
  verbose = newVerbose;
}

void VirialAlpha::setAlpha(double aCenter, double aSpan) {
  refMeter.setAlpha(aCenter, aSpan);
  refAverage.reset();
  targetMeter.setAlpha(1/aCenter, -aSpan);
  targetAverage.reset();
}

void VirialAlpha::analyze(double &jBest) {
  newAlpha = 0;
  newAlphaErr = 0;
  int numAlpha = refMeter.getNumAlpha();
  double lnRatio[numAlpha];
  const double *alpha = refMeter.getAlpha();
  double **refOverStats = refAverage.getStatistics();
  long tbc = targetAverage.getBlockCount();
  double **targetOverStats = targetAverage.getStatistics();
  for (int j=0; j<numAlpha; j++) {
    lnRatio[j] = log(refOverStats[j][AVG_AVG]/targetOverStats[j][AVG_AVG]);
    if (j>0 && lnRatio[j]*lnRatio[j-1] <= 0) {
      // linear interpolation on log scale
      double xj = lnRatio[j-1]/(lnRatio[j-1]-lnRatio[j]);
      jBest = j-1 + xj;
      newAlpha = exp(log(alpha[j-1]) + xj*(log(alpha[j]/alpha[j-1])));

      double ratio1 = exp(lnRatio[j-1]);
      double err1 = AverageRatio::ratioErr(refOverStats[j-1][AVG_AVG], refOverStats[j-1][AVG_ERR],
                               targetOverStats[j-1][AVG_AVG], targetOverStats[j-1][AVG_ERR], 0);
      double ratio2 = exp(lnRatio[j]);
      double err2 = AverageRatio::ratioErr(refOverStats[j][AVG_AVG], refOverStats[j][AVG_ERR],
                               targetOverStats[j][AVG_AVG], targetOverStats[j][AVG_ERR], 0);
      newAlphaErr = (err1/ratio1 > err2/ratio2 ? err1/ratio1 : err2/ratio2)*newAlpha;
      double ac1 = targetOverStats[j-1][AVG_ACOR];
      double ac2 = targetOverStats[j][AVG_ACOR];
      alphaCor = ac1 > ac2 ? ac1 : ac2;
      return;
    }
  }
  int jb = (fabs(lnRatio[0]) < fabs(lnRatio[numAlpha-1])) ? 0 : numAlpha-1;
  jBest = jb;
  newAlpha = alpha[jb];
  newAlphaErr = AverageRatio::ratioErr(refOverStats[jb][AVG_AVG], refOverStats[jb][AVG_ERR],
                               targetOverStats[jb][AVG_AVG], targetOverStats[jb][AVG_ERR], 0);
}

void VirialAlpha::getNewAlpha(double &na, double &nae, double &ac) {
  na = newAlpha;
  nae = newAlphaErr;
  ac = alphaCor;
}

void VirialAlpha::run() {
  while (!allDone) {
    runSteps(1000);
  }
}

void VirialAlpha::runSteps(int numSteps) {
  refIntegrator.doSteps(numSteps);
  targetIntegrator.doSteps(numSteps);
  stepCount += numSteps;
  if (stepCount >= nextCheck) {
    double jBestAlpha;
    analyze(jBestAlpha);
    if (verbose) printf("alpha  avg: %22.15e   err: %12.5e   cor: % 6.4f\n", newAlpha, newAlphaErr, alphaCor);
    int numAlpha = refMeter.getNumAlpha();
    const double* alpha = refMeter.getAlpha();
    double span = log(alpha[numAlpha-1]/alpha[0]);
    int nextCheckFac = 1.4;
    if (jBestAlpha<numAlpha*0.1 || jBestAlpha>(numAlpha-1)*0.9) span *= 2;
    else if (alphaCor < 0.3 && span > 0.5 && jBestAlpha>numAlpha*0.2 && jBestAlpha<(numAlpha-1)*0.8) span *= 0.25;
    else if (alphaCor < 0.6 && span > 0.5 && jBestAlpha>numAlpha*0.2 && jBestAlpha<(numAlpha-1)*0.8) span *= 0.6;
    else if (alphaCor > 0.2) nextCheckFac *= 2;
    else if (alphaCor < 0.1 && newAlphaErr/newAlpha < 0.02) allDone = true;
    refMeter.setAlpha(newAlpha, span);
    targetMeter.setAlpha(1/newAlpha, -span);
    refAverage.reset();
    targetAverage.reset();
    nextCheck *= nextCheckFac;
    nextCheck = stepCount + nextCheck;
  }
}
