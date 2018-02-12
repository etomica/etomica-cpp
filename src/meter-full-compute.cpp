#include "meter.h"

MeterFullCompute::MeterFullCompute(PotentialMaster& p) : Meter(0), potentialMaster(p), data(nullptr) {}

void MeterFullCompute::addCallback(PotentialCallback* pcb) {
  callbacks.push_back(pcb);
  nData += pcb->getNumData();
  data = (double*)realloc(data, nData*sizeof(double));
}

double* MeterFullCompute::getData() {
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    (*it)->reset();
  }
  potentialMaster.computeAll(callbacks);
  int idx = 0;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    double* d = (*it)->getData();
    int nd = (*it)->getNumData();
    std::copy(d, d+nd, data+idx);
    idx += nd;
  }
  return data;
}

