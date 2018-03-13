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
    long getBlockSize() {return blockSize;}
    long getBlockCount() {return blockCount;}
    void setNumData(int newNumData);
    virtual void reset();
};

// also computes ratio of each quantity with the last quantity
class AverageRatio : public Average {
  private:
    static double ratioErr(double nAvg, double nErr, double dAvg, double dErr, double cor) {
      if (nAvg==0 && nErr==0) return 0;
      double ratio = nAvg/dAvg;
      if (nAvg==0) {
        return sqrt((nErr*nErr)/(dAvg*dAvg));
      }
      return sqrt((nErr*nErr/(nAvg*nAvg) + dErr*dErr/(dAvg*dAvg) - 2*cor*nErr*dErr/(nAvg*dAvg)) * ratio*ratio);
    }

  protected:
    double **ratioStats;
    double **ratioCovariance;
  public:
    AverageRatio(int nData, long blockSize, long maxBlockCount, bool doCovariance);
    ~AverageRatio();
    virtual void reset();
    double** getRatioStatistics();
    double** getRatioCovariance();
};
