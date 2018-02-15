#include "average.h"

DataSink::DataSink() : pushInterval(1), intervalCountdown(1) {}

void DataSink::setPushInterval(int interval) {
  pushInterval = interval;
  intervalCountdown = interval;
}

void DataSink::addDataSink(DataSink* s) {
  sinks.push_back(s);
}


