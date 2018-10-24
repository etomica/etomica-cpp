/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "cluster.h"

Cluster::Cluster(int nm, int nv, bool cached) : IntegratorListener(), nValues(nv), numMolecules(nm), useCache(cached), cacheDirty(true), inTrial(false) {
  callAccept = true;
  callReject = true;
  values = new double[nValues];
  oldValues = new double[nValues];
}

Cluster::~Cluster() {
  delete[] values;
  delete[] oldValues;
}

void Cluster::moveRejected(MCMove& move, double chi) {
  if (!useCache || !inTrial) return;
  trialRejected();
}

void Cluster::moveAccepted(MCMove& move, double chi) {
  if (!useCache) return;
  if (!inTrial) {
    cacheDirty = true;
    return;
  }
  inTrial = false;
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
  inTrial = true;
}

void Cluster::trialRejected() {
  for (int m=0; m<nValues; m++) values[m] = oldValues[m];
  cacheDirty = false;
  inTrial = false;
}

