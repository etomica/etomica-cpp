/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "integrator.h"
#include "meter.h"

class DataSink;

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
