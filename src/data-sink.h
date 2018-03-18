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

  public:
    Average(int nData, long blockSize, long maxBlockCount, bool doCovariance);
    virtual ~Average();
    virtual void addData(double* x);
    double** getStatistics();
    double** getBlockCovariance();
    double** getBlockCorrelation();
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
      return ((vi/vd)/eid)*((vj/vd)/ejd)*((ei/vi)*(ed/vd) + (ei/vi)*(ej/vj)*cij - (ei/vi)*(ed/vd)*cid - (ej/vj)*(ed/vj)*cjd);
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
      return ((vi/vd)/eid)*((vj/vd)/ejd)*((ei/vi)*(ed/vd) + (ei/vi)*(ej/vj)*cij - (ei/vi)*(ed/vd)*cid - (ej/vj)*(ed/vj)*cjd)/(eid*ejd);
    }
};
