#include "meter.h"

DataPump::DataPump(Meter& m, int i) : IntegratorListener(), meter(m), interval(i), intervalCountdown(i) {
  sink = new Average(m.getNumData(),10,100);
}

DataPump::DataPump(Meter& m, int i, DataSink* s) : IntegratorListener(), meter(m), interval(i), intervalCountdown(i), sink(s) {}

DataPump::~DataPump() {}

void DataPump::stepFinished() {
  if (--intervalCountdown > 0) return;
  sink->addData(meter.getData());
  intervalCountdown = interval;
}

DataSink* DataPump::getDataSink() {
  return sink;
}
