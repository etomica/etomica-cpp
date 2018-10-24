/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
    double icor = d <= 0 ? 0 : blockCovariance[i][nData-1] / sqrt(d);
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
  double vd = stats[nData-1][AVG_AVG];
  double ed = stats[nData-1][AVG_ERR];
  for (int i=0; i<nData; i++) {
    double vi = stats[i][AVG_AVG];
    double ei = stats[i][AVG_ERR];
    double x = blockCovariance[i][i] * blockCovariance[nData-1][nData-1];
    if (x <= 0) {
      for (int j=0; j<=i; j++) ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
      continue;
    }
    double cid = blockCovariance[i][nData-1] / sqrt(x);
    for (int j=0; j<=i; j++) {
      double vj = stats[j][AVG_AVG];
      double ej = stats[j][AVG_ERR];
      if (blockCovariance[j][j] <= 0) {
        ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
        continue;
      }
      double cjd = blockCovariance[j][nData-1] / sqrt(blockCovariance[j][j] * blockCovariance[nData-1][nData-1]);
      double cij = blockCovariance[i][j] / sqrt(blockCovariance[i][i]*blockCovariance[j][j]);
      ratioCovariance[j][i] = ratioCovariance[i][j] = ratioCov(vi, vj, vd, ei, ej, ed, cij, cid, cjd);
    }
  }
  return ratioCovariance;
}

double** AverageRatio::getRatioCorrelation() {
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
  double vd = stats[nData-1][AVG_AVG];
  double ed = stats[nData-1][AVG_ERR];
  for (int i=0; i<nData; i++) {
    double vi = stats[i][AVG_AVG];
    double ei = stats[i][AVG_ERR];
    double x = blockCovariance[i][i] * blockCovariance[nData-1][nData-1];
    if (x <= 0) {
      for (int j=0; j<=i; j++) ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
      continue;
    }
    double cid = blockCovariance[i][nData-1] / sqrt(x);
    for (int j=0; j<=i; j++) {
      double vj = stats[j][AVG_AVG];
      double ej = stats[j][AVG_ERR];
      if (blockCovariance[j][j] <= 0) {
        ratioCovariance[j][i] = ratioCovariance[i][j] = 0;
        continue;
      }
      double cjd = blockCovariance[j][nData-1] / sqrt(blockCovariance[j][j] * blockCovariance[nData-1][nData-1]);
      double cij = blockCovariance[i][j] / sqrt(blockCovariance[i][i]*blockCovariance[j][j]);
      ratioCovariance[j][i] = ratioCovariance[i][j] = ratioCor(vi, vj, vd, ei, ej, ed, cij, cid, cjd);
    }
  }
  return ratioCovariance;
}
