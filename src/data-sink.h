#pragma once

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
  private:
    int nData;
    long defaultBlockSize, blockSize, blockCount, maxBlockCount;
    long blockCountdown;
    double *mostRecent;
    double *currentBlockSum, *blockSum, *blockSum2, *correlationSum;
    double** stats;
    double** blockSums;
    double** blockCovariance;
    double** blockCovSum;
    const bool doCovariance;

    void collapseBlocks();

  public:
    Average(int nData, long blockSize, long maxBlockCount, bool doCovariance);
    ~Average();
    void addData(double* x);
    double** getStatistics();
    double** getBlockCovariance();
    long getBlockSize() {return blockSize;}
    long getBlockCount() {return blockCount;}
    void setNumData(int newNumData);
    void reset();
};
