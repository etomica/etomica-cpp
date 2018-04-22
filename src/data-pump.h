#include "integrator.h"
#include "meter.h"

class DataPump : public IntegratorListener {
  private:
    Meter& meter;
    int interval;
    int intervalCountdown;
    vector<DataSink*> sinks;
    const bool dataSink1isMine;
  public:
    DataPump(Meter& meter, int interval);
    DataPump(Meter& meter, int interval, DataSink* sink);
    ~DataPump();
    void stepFinished();
    DataSink* getDataSink(int i);
    void addDataSink(DataSink* sink);
};
