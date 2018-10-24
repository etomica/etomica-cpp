/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "meter.h"
#include "cluster.h"

class MeterVirialDirect : public Meter {
  private:
    Cluster &targetCluster;
    Cluster &refCluster;
    double *data;
  public:
    MeterVirialDirect(Cluster &targetCluster, Cluster &refCluster);
    ~MeterVirialDirect();
    double* getData();
};

class MeterVirialOverlap : public Meter {
  private:
    Cluster &primaryCluster;
    Cluster &perturbCluster;
    double *data;
    double *alpha;
    const int numAlpha;

  public:
    MeterVirialOverlap(Cluster &clusterPrimary, Cluster &clusterPerturb, double alphaCenter, double alphaSpan, int numAlpha);
    ~MeterVirialOverlap();
    void setAlpha(double alphaCenter, double alphaSpan);
    int getNumAlpha();
    const double* getAlpha();
    double* getData();
};
