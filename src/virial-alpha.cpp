/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "virial.h"

// perhaps just take Clusters and make Meter and Average internally

VirialAlpha::VirialAlpha(IntegratorMC &rIntegrator, IntegratorMC &tIntegrator, Cluster &refClusterRef, Cluster &refClusterTarget, Cluster &targetClusterRef, Cluster &targetClusterTarget) :
    stepCount(0), nextCheck(1000), refIntegrator(rIntegrator), targetIntegrator(tIntegrator),
    refMeter(MeterVirialOverlap(refClusterRef, refClusterTarget, 1, 5, 11)),
    targetMeter(MeterVirialOverlap(targetClusterTarget, targetClusterRef, 1, -5, 11)),
    refAverage(11, 1, 1000, false), targetAverage(11, 1, 100, false),
    refPump(refMeter,1,&refAverage), targetPump(targetMeter,1,&targetAverage),
    newAlpha(0), newAlphaErr(0), alphaCor(0), alphaSpan(0), allDone(false), verbose(false), disposed(false) {
  refIntegrator.addListener(&refPump);
  targetIntegrator.addListener(&targetPump);
  int numAlpha = refMeter.getNumAlpha();
  const double *alpha = refMeter.getAlpha();
  alphaSpan = log(alpha[numAlpha-1]/alpha[0]);
}

VirialAlpha::~VirialAlpha() {
  dispose();
}

void VirialAlpha::dispose() {
  if (disposed) return;
  refIntegrator.removeListener(&refPump);
  targetIntegrator.removeListener(&targetPump);
  for (vector<double*>::iterator it = allAlphaStats.begin(); it!=allAlphaStats.end(); it++) {
    free(*it);
  }
  allAlphaStats.clear();
  disposed = true;
}

int VirialAlpha::getNumSavedStats() {
  return allAlphaStats.size();
}

double* VirialAlpha::getSavedStats(int i) {
  return allAlphaStats[i];
}

void VirialAlpha::setVerbose(bool newVerbose) {
  verbose = newVerbose;
}

void VirialAlpha::setAlpha(double aCenter, double aSpan) {
  alphaSpan = aSpan;
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

double* VirialAlpha::getAlphaStatistics() {
  alphaStats[0] = newAlpha;
  alphaStats[1] = newAlphaErr;
  alphaStats[2] = alphaCor;
  alphaStats[3] = alphaSpan;
  return alphaStats;
}

void VirialAlpha::run() {
  while (!allDone) {
    runSteps(1000);
  }
}

bool VirialAlpha::getAllDone() {
  return allDone;
}

void VirialAlpha::runSteps(int numSteps) {
  refIntegrator.doSteps(numSteps);
  targetIntegrator.doSteps(numSteps);
  stepCount += numSteps;
  if (stepCount >= nextCheck) {
    double jBestAlpha;
    analyze(jBestAlpha);
    double* saveStats = (double*)malloc(6*sizeof(double));
    // both ref and target avg have the same counts
    saveStats[0] = refAverage.getBlockCount() * refAverage.getBlockSize();
    saveStats[1] = refAverage.getBlockSize();
    saveStats[2] = newAlpha;
    saveStats[3] = newAlphaErr;
    saveStats[4] = alphaCor;
    saveStats[5] = alphaSpan;
    allAlphaStats.push_back(saveStats);
    if (verbose) printf("alpha  avg: %22.15e   err: %12.5e   cor: % 6.4f\n", newAlpha, newAlphaErr, alphaCor);
    int numAlpha = refMeter.getNumAlpha();
    double nextCheckFac = 1.4;
    if (jBestAlpha<numAlpha*0.1 || jBestAlpha>(numAlpha-1)*0.9) alphaSpan *= 2;
    else if (alphaCor < 0.3 && alphaSpan > 0.5 && jBestAlpha>numAlpha*0.2 && jBestAlpha<(numAlpha-1)*0.8) alphaSpan *= 0.25;
    else if (alphaCor < 0.6 && alphaSpan > 0.5 && jBestAlpha>numAlpha*0.2 && jBestAlpha<(numAlpha-1)*0.8) alphaSpan *= 0.6;
    else if (alphaCor > 0.2) nextCheckFac *= 2;
    else if (alphaCor < 0.1 && newAlphaErr/newAlpha < 0.02) allDone = true;
    setAlpha(newAlpha, alphaSpan);
    nextCheck *= nextCheckFac;
    nextCheck = stepCount + nextCheck;
  }
}
