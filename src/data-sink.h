/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <math.h>
#include <vector>

using namespace std;

class DataSink {
  protected:
    vector<DataSink*> sinks;
    int pushInterval, intervalCountdown;
  public:
    DataSink();
    virtual ~DataSink() {}
    virtual void addData(double* x) = 0;
    void addDataSink(DataSink* sink);
    void setPushInterval(int newPushInterval);
};

#define AVG_CUR 0
#define AVG_AVG 1
#define AVG_ERR 2
#define AVG_ACOR 3

class Average : public DataSink {
  protected:
    int nData;
    long defaultBlockSize, blockSize, blockCount, maxBlockCount;
    long blockCountdown;
    double *mostRecent;
    double *currentBlockSum, *blockSum, *blockSum2, *correlationSum;
    double *prevBlockSum, *firstBlockSum;
    double **stats;
    double **blockSums;
    double **blockCovariance;
    double **blockCovSum;
    const bool doCovariance;

    void collapseBlocks();

    void dispose();
    void unset();

  public:
    Average(int nData, long blockSize, long maxBlockCount, bool doCovariance);
    virtual ~Average();
    virtual void addData(double* x);
    double** getStatistics();
    double** getBlockCovariance();
    double** getBlockCorrelation();
    void setBlockSize(long blockSize);
    void setMaxBlockCount(long maxBlockCount);
    long getBlockSize() {return blockSize;}
    long getBlockCount() {return blockCount;}
    void setNumData(int newNumData);
    int getNumData();
    virtual void reset();
};

// also computes ratio of each quantity with the last quantity
class AverageRatio : public Average {

  protected:
    double **ratioStats;
    double **ratioCovariance;
  public:
    AverageRatio(int nData, long blockSize, long maxBlockCount, bool doCovariance);
    ~AverageRatio();
    virtual void reset();
    double** getRatioStatistics();
    double** getRatioCovariance();
    double** getRatioCorrelation();

    static double ratioErr(double nAvg, double nErr, double dAvg, double dErr, double cor) {
      if (nAvg==0 && nErr==0) return 0;
      double ratio = nAvg/dAvg;
      if (nAvg==0) {
        return sqrt((nErr*nErr)/(dAvg*dAvg));
      }
      return sqrt((nErr*nErr/(nAvg*nAvg) + dErr*dErr/(dAvg*dAvg) - 2*cor*nErr*dErr/(nAvg*dAvg)) * ratio*ratio);
    }

    /**
     * Compute covariance of i/d and j/d
     *
     * vi: value of i
     * vj: value of j
     * vd: value of d
     * ei: error in the i numerator
     * ej: error in the j numerator
     * ed: error in the denominator d
     * cij: correlation between i and j
     * cid: correlation between i and d
     * cjd: correlation between j and d
     */
    static double ratioCov(double vi, double vj, double vd, double ei, double ej, double ed, double cij, double cid, double cjd) {
      double eid = ratioErr(vi, ei, vd, ed, cid);
      double ejd = ratioErr(vj, ej, vd, ed, cjd);
      if (eid==0 || ejd==0) return 0;
      return (vi/vd)*(vj/vd)*((ei/vi)*(ed/vd) + (ei/vi)*(ej/vj)*cij - (ei/vi)*(ed/vd)*cid - (ej/vj)*(ed/vj)*cjd);
    }

    /**
     * Compute correlation of i/d and j/d
     *
     * vi: value of i
     * vj: value of j
     * vd: value of d
     * ei: error in the i numerator
     * ej: error in the j numerator
     * ed: error in the denominator d
     * cij: correlation between i and j
     * cid: correlation between i and d
     * cjd: correlation between j and d
     */
    static double ratioCor(double vi, double vj, double vd, double ei, double ej, double ed, double cij, double cid, double cjd) {
      double eid = ratioErr(vi, ei, vd, ed, cid);
      double ejd = ratioErr(vj, ej, vd, ed, cjd);
      if (eid==0 || ejd==0) return 0;
      return ((vi/vd)/eid)*((vj/vd)/ejd)*((ed/vd)*(ed/vd) + (ei/vi)*(ej/vj)*cij - (ei/vi)*(ed/vd)*cid - (ej/vj)*(ed/vd)*cjd);
    }
};

class History : public DataSink {
  protected:
    int nData, historySize;
    // 0: scrolling, 1: complete, 2: collapsing discard, 3: collapsing avg
    int historyType;
    int count, pointer, collapseSize, skipCount;
    double *collapseSum;
    double **data, **history;

    void collapseDiscard();
    void collapseAverage();
    void dispose();
    void unset();

  public:
    History(int nData, int historyType, int historySize);
    virtual ~History();
    virtual void addData(double* x);
    double** getHistory();
    long getRawCount() {return count;}
    int getHistorySize();
    void setNumData(int newNumData);
    int getNumData();
    virtual void reset();
};

