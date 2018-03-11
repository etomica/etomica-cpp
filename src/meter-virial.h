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
