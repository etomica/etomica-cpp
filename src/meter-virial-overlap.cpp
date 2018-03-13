#include <cmath>
#include "meter-virial.h"

// This meter takes 2 clusters -- a primary and a perturb cluster.  The sampling weight (pi) is
// taken to be the absolute value of the first value from the primary cluster.  The absolute
// value of the first value from the perturb cluster is used as the sampling weight (piP) of the 
// other system (reference or target).  The overlap average is formed from
//
// gamma_OS = piT piR / (piT + alpha piR)
//
// where T(arget) or R(eference) could be the primary or perturb cluster.  The meter computes
// |gamma_OS / pi|.  When the primary cluster is the reference cluster, then the alpha given to
// the meter should be 1/alpha and the alpha span should be negative.  The meter computes
// |alpha gamma_OS / pi| for this case.
//
// This meter can operate with multiple alpha values or a single alpha value.  With multiple
// alpha values (used when trying to find the best alpha), the meter measures only the overlap
// averages with the various alpha values.  With a single alpha value (used when running the
// production stage of the calculation), the first n values are v_i/pi, where v_i is the ith
// value from the primary cluster, which returns n values.  The last value returned by the
// meter is the overlap ratio as described above.

MeterVirialOverlap::MeterVirialOverlap(Cluster &cluster1, Cluster &cluster2, double aCenter, double aSpan, int nAlpha) : Meter(0), primaryCluster(cluster1), perturbCluster(cluster2), numAlpha(nAlpha) {
  if (nAlpha>1) {
    nData = nAlpha;
    if (aSpan == 0) {
      fprintf(stderr, "If # of alpha > 1, then alpha span can't be 0");
      abort();
    }
  }
  else {
    nData = primaryCluster.numValues() + 1;
    if (aSpan != 0) {
      fprintf(stderr, "If # of alpha is 1, then alpha span must be 0");
      abort();
    }
  }
  data = new double[nData];
  alpha = new double[nData];
  setAlpha(aCenter, aSpan);
}

MeterVirialOverlap::~MeterVirialOverlap() {
  delete[] data;
  delete[] alpha;
}

void MeterVirialOverlap::setAlpha(double aCenter, double aSpan) {
  if (numAlpha == 1) {
    alpha[0] = aCenter;
    return;
  }
  for (int i=0; i<numAlpha; i++) {
    alpha[i] = exp((i-(numAlpha-1)/2) * aSpan);
  }
}

double* MeterVirialOverlap::getData() {
  const double* primaryValues = primaryCluster.getValues();
  double pi = fabs(primaryValues[0]);
  if (pi == 0 || pi == std::numeric_limits<double>::infinity() || std::isnan(pi)) {
    fprintf(stderr, "pi is %f\n", pi);
    abort();
  }
  double perturbValue = fabs(perturbCluster.getValues()[0]);
  if (numAlpha == 1) {
    for (int i=0; i<primaryCluster.numValues(); i++) {
      data[i] = primaryValues[i] / pi;
    }
    data[nData-1] = perturbCluster.getValues()[0] / pi;
  }
  else {
    for (int i=0; i<numAlpha; i++) {
      // gamma_OS = pi1 pi0 / (pi1 + alpha pi0)
      // 1: gamma_OS/pi1 = pi0 / (pi1 + alpha pi0)
      // 0: gamma_OS/pi0 = pi1 / (pi1 + alpha pi0)
      //                 = (1/alpha) pi1 / (pi0 + (1/alpha) pi1)
      // for 0 case, we use negative alphaSpan (alpha => 1/alpha) and compute: alpha gammaOS/pi0
      data[i] = perturbValue / (pi + alpha[i]*perturbValue);
    }
  }
  return data;
}
