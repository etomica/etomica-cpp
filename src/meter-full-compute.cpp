/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "meter.h"
#include "potential-master.h"

MeterFullCompute::MeterFullCompute(PotentialMaster& p) : Meter(0), potentialMaster(p), data(nullptr), doCompute(true) {}

MeterFullCompute::~MeterFullCompute() {
  free(data);
}

void MeterFullCompute::addCallback(PotentialCallback* pcb) {
  callbacks.push_back(pcb);
  nData += pcb->getNumData();
  data = (double*)realloc(data, nData*sizeof(double));
}

void MeterFullCompute::setDoCompute(bool doC) {
  doCompute = doC;
}

double* MeterFullCompute::getData() {
  if (doCompute) {
    // probably MC, we need to invoke computeAll ourselves
    // MD can invoke callbacks internally every step, then we just retrieve the data
    for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
      (*it)->reset();
    }
    potentialMaster.computeAll(callbacks);
  }
  int idx = 0;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    double* d = (*it)->getData();
    int nd = (*it)->getNumData();
    std::copy(d, d+nd, data+idx);
    idx += nd;
  }
  return data;
}

