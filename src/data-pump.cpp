#include "data-pump.h"

DataPump::DataPump(Meter& m, int i) : IntegratorListener(), meter(m), interval(i), intervalCountdown(i), dataSink1isMine(true) {
  callStepFinished = true;
  DataSink *s = new Average(m.getNumData(),1,200,true);
  addDataSink(s);
}

DataPump::DataPump(Meter& m, int i, DataSink* s) : IntegratorListener(), meter(m), interval(i), intervalCountdown(i), dataSink1isMine(false) {
  callStepFinished = true;
  addDataSink(s);
}

DataPump::~DataPump() {
  if (dataSink1isMine) {
    delete sinks[0];
  }
}

void DataPump::stepFinished() {
  if (--intervalCountdown > 0) return;
  double* mData = meter.getData();
  for (vector<DataSink*>::iterator it = sinks.begin(); it!=sinks.end(); it++) {
    (*it)->addData(mData);
  }
  intervalCountdown = interval;
}

DataSink* DataPump::getDataSink(int i) {
  return sinks[i];
}

void DataPump::addDataSink(DataSink* s) {
  sinks.push_back(s);
}
