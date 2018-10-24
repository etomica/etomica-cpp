/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <cmath>
#include "meter-virial.h"

MeterVirialDirect::MeterVirialDirect(Cluster &tCluster, Cluster &rCluster) : Meter(0), targetCluster(tCluster), refCluster(rCluster) {
  nData = targetCluster.numValues() + 1;
  data = new double[nData];
}

MeterVirialDirect::~MeterVirialDirect() {
  delete[] data;
}

double* MeterVirialDirect::getData() {
  const double* tarValues = targetCluster.getValues();
  double pi = fabs(tarValues[0]);
  if (pi == 0 || pi == std::numeric_limits<double>::infinity() || std::isnan(pi)) {
    fprintf(stderr, "pi is %f\n", pi);
    abort();
  }
  for (int j=0; j<nData-1; j++) {
    data[j] = tarValues[j] / pi;
  }
  data[nData-1] = refCluster.getValues()[0] / pi;
  return data;
}
