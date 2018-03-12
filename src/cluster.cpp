#include "cluster.h"

Cluster::Cluster(int nm, int nv, bool cached) : nValues(nv), numMolecules(nm), useCache(cached), cacheDirty(true) {
  values = new double[nValues];
  oldValues = new double[nValues];
}

Cluster::~Cluster() {
  delete[] values;
  delete[] oldValues;
}

void Cluster::setCachingEnabled(bool enabled) {
  useCache = enabled;
}

int Cluster::numValues() {
  return nValues;
}

void Cluster::trialNotify() {
  for (int m=0; m<nValues; m++) oldValues[m] = values[m];
  cacheDirty = true;
}

void Cluster::trialRejected() {
  for (int m=0; m<=nValues; m++) values[m] = oldValues[m];
  cacheDirty = false;
}

