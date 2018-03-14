#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "data-sink.h"
#include "alloc2d.h"

Average::Average(int n, long bs, long mBC, bool doCov) : nData(n), defaultBlockSize(bs), blockSize(bs), blockCount(0),
                               maxBlockCount(mBC), blockCountdown(bs), doCovariance(doCov) {
  mostRecent = nullptr;
  currentBlockSum = nullptr;
  blockSum = nullptr;
  blockSum2 = nullptr;
  correlationSum = nullptr;
  blockSums = nullptr;
  stats = nullptr;
  blockCovariance = nullptr;
  blockCovSum = nullptr;
  reset();
}

Average::~Average() {
  free(mostRecent);
  free(currentBlockSum);
  free(blockSum);
  free(blockSum2);
  free(correlationSum);
  free2D((void**)blockSums);
  free2D((void**)stats);
  if (doCovariance) {
    free2D((void**)blockCovariance);
    free2D((void**)blockCovSum);
  }
}

void Average::setNumData(int newNumData) {
  nData = newNumData;
  reset();
}

int Average::getNumData() {
  return nData;
}

void Average::reset() {
  blockCount = 0;
  blockSize = defaultBlockSize;
  blockCountdown = blockSize;
  // realloc our arrays so that we can adjust if n changes
  mostRecent = (double*)realloc(mostRecent, nData*sizeof(double));
  currentBlockSum = (double*)realloc(currentBlockSum, nData*sizeof(double));
  blockSum = (double*)realloc(blockSum, nData*sizeof(double));
  blockSum2 = (double*)realloc(blockSum2, nData*sizeof(double));
  correlationSum = (double*)realloc(correlationSum, nData*sizeof(double));
  for (int i=0; i<nData; i++) {
    currentBlockSum[i] = blockSum[i] = blockSum2[i] = correlationSum[i] = 0;
  }
  if (maxBlockCount>0) {
    if (maxBlockCount%2==1 || maxBlockCount < 4) {
      fprintf(stderr, "Not nice!  Give me a max block count that's even and >= 4!\n");
      exit(0);
    }
    blockSums = (double**)realloc2D((void**)blockSums, nData, maxBlockCount, sizeof(double));
  }
  stats = (double**)realloc2D((void**)stats, nData, 4, sizeof(double));
  if (doCovariance) {
    blockCovSum = (double**)realloc2D((void**)blockCovSum, nData, nData, sizeof(double));
    for (int i=0; i<nData; i++) {
      for (int j=0; j<=i; j++) {
        blockCovSum[i][j] = 0;
      }
    }
    blockCovariance = (double**)realloc2D((void**)blockCovariance, nData, nData, sizeof(double));
  }
}

void Average::addData(double *x) {
  for (int i=0; i<nData; i++) {
    mostRecent[i] = x[i];
    currentBlockSum[i] += x[i];
  }
  if (--blockCountdown == 0) {
    if (doCovariance) {
      double blockSizeSq = blockSize*blockSize;
      for (int i=0; i<nData; i++) {
        for (int j=0; j<=i; j++) {
          double ijx = currentBlockSum[i]*currentBlockSum[j]/blockSizeSq;
          blockCovSum[i][j] += ijx;
        }
      }
    }
    for (int i=0; i<nData; i++) {
      blockSum[i] += currentBlockSum[i];
      currentBlockSum[i] /= blockSize;
      if (blockCount>0) correlationSum[i] += blockSums[i][blockCount-1] * currentBlockSum[i];
      blockSums[i][blockCount] = currentBlockSum[i];
      blockSum2[i] += currentBlockSum[i] * currentBlockSum[i];
      currentBlockSum[i] = 0;
    }
    blockCount++;
    blockCountdown = blockSize;

    if (blockCount == maxBlockCount) {
      collapseBlocks();
    }
  }
  // we don't push our data
}

void Average::collapseBlocks() {
  blockCount /= 2;
  for (int i=0; i<nData; i++) {
    blockSum2[i] = 0;
    correlationSum[i] = 0;
    // the first half of the blocks contain all previous data
    if (doCovariance) {
      for (int k=0; k<=i; k++) {
        blockCovSum[i][k] = 0;
      }
    }
    for (int j=0; j<blockCount; j++) {
      blockSums[i][j] = (blockSums[i][2*j] + blockSums[i][2*j+1]) / 2;
      blockSum2[i] += blockSums[i][j] * blockSums[i][j];
      if (j>0) {
        correlationSum[i] += blockSums[i][j-1] * blockSums[i][j];
      }
      if (doCovariance) {
        for (int k=0; k<=i; k++) {
          blockCovSum[i][k] += blockSums[i][j]*blockSums[k][j];
        }
      }
    }
    for (int j=blockCount; j<maxBlockCount; j++) {
      blockSums[i][j] = 0;
    }
  }
  blockSize *= 2;
  blockCountdown = blockSize;
  getBlockCovariance();
}

double** Average::getStatistics() {
  if (blockCount==0) {
    for (int i=0; i<nData; i++) {
      stats[i][AVG_CUR] = mostRecent[i];
      stats[i][1] = stats[i][2] = stats[i][3] = NAN;
    }
    return stats;
  }
  for (int i=0; i<nData; i++) {
    stats[i][AVG_CUR] = mostRecent[i];
    stats[i][AVG_AVG] = blockSum[i] / (blockSize*blockCount);
    if (blockCount == 1) {
      for (int i=0; i<nData; i++) {
        stats[i][AVG_ERR] = stats[i][AVG_ACOR] = NAN;
      }
      continue;
    }
    stats[i][AVG_ERR] = blockSum2[i] / blockCount - stats[i][AVG_AVG]*stats[i][AVG_AVG];
    if (stats[i][AVG_ERR]<0) stats[i][AVG_ERR] = 0;
    if (stats[i][AVG_ERR] == 0) {
      stats[i][AVG_ACOR] = 0;
    }
    else {
      double bc = (((2 * blockSum[i] / blockSize - blockSums[i][0] - blockSums[i][blockCount-1]) * stats[i][AVG_AVG] - correlationSum[i]) / (1-blockCount) + stats[i][AVG_AVG]*stats[i][AVG_AVG]) / stats[i][AVG_ERR];
      stats[i][AVG_ACOR] = (isnan(bc) || bc <= -1 || bc >= 1) ? 0 : bc;
    }
    stats[i][AVG_ERR] = sqrt(stats[i][AVG_ERR]/(blockCount-1));
  }
  return stats;
}

double** Average::getBlockCovariance() {
  double totSamples = blockSize*blockCount;
  double totSq = totSamples*totSamples;
  for (int i=0; i<nData; i++) {
    for (int j=0; j<=i; j++) {
      blockCovariance[i][j] = blockCovariance[j][i] = blockCovSum[i][j]/blockCount - blockSum[i]*blockSum[j]/totSq;
    }
  }

  return blockCovariance;
}
