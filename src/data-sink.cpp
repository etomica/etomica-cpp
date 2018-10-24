/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "data-sink.h"

DataSink::DataSink() : pushInterval(1), intervalCountdown(1) {}

void DataSink::setPushInterval(int interval) {
  pushInterval = interval;
  intervalCountdown = interval;
}

void DataSink::addDataSink(DataSink* s) {
  sinks.push_back(s);
}


