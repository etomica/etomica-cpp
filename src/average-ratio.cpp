#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "data-sink.h"
#include "alloc2d.h"

AverageRatio::AverageRatio(int n, long bs, long mBC, bool doCov) : Average(n, bs, mBC, doCov), ratioStats(nullptr), ratioCovariance(nullptr) {
  reset();
}

AverageRatio::~AverageRatio() {
  free2D((void**)ratioStats);
  if (doCovariance) free2D((void**)ratioCovariance);
}

void AverageRatio::reset() {
  Average::reset();
  ratioStats = (double**)realloc2D((void**)ratioStats, nData, 3, sizeof(double));
  if (doCovariance) ratioCovariance = (double**)realloc2D((void**)ratioCovariance, nData, nData, sizeof(double));
}

double** AverageRatio::getRatioStatistics() {
  if (blockCount==0) {
    for (int i=0; i<nData; i++) {
      ratioStats[i][AVG_CUR] = NAN;
      ratioStats[i][AVG_AVG] = ratioStats[i][AVG_ERR] = NAN;
    }
    return ratioStats;
  }
  getStatistics();
  getBlockCovariance();
  for (int i=0; i<nData; i++) {
    ratioStats[i][AVG_CUR] = mostRecent[i]/mostRecent[nData-1];
    ratioStats[i][AVG_AVG] = blockSum[i] / blockSum[nData-1];
    if (blockCount == 1) {
      for (int i=0; i<nData; i++) {
        ratioStats[i][AVG_ERR] = NAN;
      }
      continue;
    }
    double d = blockCovariance[i][i] * blockCovariance[nData-1][nData-1];
    double icor = d <= 0 ? 0 : blockCovariance[i][nData-1] / sqrt(blockCovariance[i][i] * blockCovariance[nData-1][nData-1]);
    ratioStats[i][AVG_ERR] = ratioErr(stats[i][AVG_AVG], stats[i][AVG_ERR], stats[nData-1][AVG_AVG], stats[nData-1][AVG_ERR], icor);
  }
  return ratioStats;
}

double** AverageRatio::getRatioCovariance() {
  if (blockCount<2) {
    for (int i=0; i<nData; i++) {
      for (int j=0; j<nData; j++) {
        ratioCovariance[i][j] = NAN;
      }
    }
    return ratioCovariance;
  }
  getStatistics();
  getBlockCovariance();
  double d = ratioStats[nData-1][AVG_ERR] / ratioStats[nData-1][AVG_AVG];
  for (int i=0; i<nData; i++) {
    double ei = ratioStats[i][AVG_ERR] / ratioStats[i][AVG_AVG];
    double ci = blockCovariance[i][nData-1] / sqrt(blockCovariance[i][i] * blockCovariance[nData-1][nData-1]);
    for (int j=0; j<=i; j++) {
      double ej = ratioStats[j][AVG_ERR] / ratioStats[j][AVG_AVG];
      double cj = blockCovariance[j][nData-1] / sqrt(blockCovariance[j][j] * blockCovariance[nData-1][nData-1]);
      double cij = blockCovariance[i][j] / sqrt(blockCovariance[i][i]*blockCovariance[j][j]);
      ratioCovariance[j][i] = ratioCovariance[i][j] = d*d + ei*ej*cij - ei*d*ci - ej*d*cj;
    }
  }
  return stats;
}
