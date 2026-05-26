/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <cmath>
#include <limits>
#include "meter-virial.h"
#include "vector.h"

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
  for (int i=0; i<100; i++)
  {
    histogram[i] = 0;
  }
  counter = 0;
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
    alpha[i] = aCenter * exp((i-(numAlpha-1.0)/2) / (numAlpha-1.0) * aSpan);
  }
}

int MeterVirialOverlap::getNumAlpha() {
  return numAlpha;
}

const double* MeterVirialOverlap::getAlpha() {
  return alpha;
}

double* MeterVirialOverlap::getData()
{
  const double* primaryValues = primaryCluster.getValues();
  double pi = fabs(primaryValues[0]);
  if (pi == 0 || pi == std::numeric_limits<double>::infinity() || std::isnan(pi)) {
    fprintf(stderr, "pi is %f in %d\n", pi, idx);
    abort();
  }
  double perturbValue = fabs(perturbCluster.getValues()[0]);
  if (numAlpha == 1) {
    for (int i=0; i<primaryCluster.numValues(); i++) {
      data[i] = primaryValues[i] / pi;
    }
    data[nData-1] = perturbValue / (perturbValue + alpha[0]*pi);
    // printf("%f, %f, %f \n",perturbValue, alpha[0], pi);
  }
  else {
    for (int i=0; i<numAlpha; i++) {
      // gamma_OS = pi1 pi0 / (pi1 + alpha pi0)
      // 0: gamma_OS/pi0 = pi1 / (pi1 + alpha pi0)
      // 1: gamma_OS/pi1 = pi0 / (pi1 + alpha pi0)
      //                 = (1/alpha) pi0 / (pi1 + (1/alpha) pi0)
      // for 1 case, we use negative alphaSpan (alpha => 1/alpha)
      //    and effectively compute: alpha gammaOS/pi1
      // <0>/<1> = (1/alpha) <gammaOS/pi0>0 / <gammaOS/pi1>1
      //         ~= 1 (when alpha is optimal)
      data[i] = perturbValue / (perturbValue + alpha[i]*pi);
    }
  }
  // if (idx == 1)
  // {
  //   printf("%f \n", primaryValues[0]);
  // }
  static int count = 0;

  if (box)
  {
    double dr12[3], dr13[3], dr23[3], shift[3];
    for (int i=0; i<3; i++) {dr12[i] = box->getAtomPosition(0)[i] - box->getAtomPosition(3)[i];}
    for (int i=0; i<3; i++) {dr13[i] = box->getAtomPosition(0)[i] - box->getAtomPosition(6)[i];}
    for (int i=0; i<3; i++) {dr23[i] = box->getAtomPosition(3)[i] - box->getAtomPosition(6)[i];}
    double dist12 = sqrt(Vector::squared(dr12));
    double dist13 = sqrt(Vector::squared(dr13));
    double dist23 = sqrt(Vector::squared(dr23));
    double dist = std::max({dist12, dist13, dist23});
    count++;
    if (count % 500 == 0) printf("%d unscaled 12: %f 13: %f 23: %f max: %f \n", count, dist12, dist13, dist23, dist);
    for (int i=0; i<3; i++) shift[i] = 0.01 * dr12[i];
    if (count == 2999)
    {
      printf("%d unscaled 12: %f 13: %f 23: %f max: %f \n", count, dist12, dist13, dist23, dist);

      for (int m = 0; m < 100; m++)
      {
        for (int k = 3; k < 6; k++) for (int i=0; i<3; i++) box->getAtomPosition(k)[i] += shift[i];

        for (int i=0; i<3; i++) {dr12[i] = box->getAtomPosition(0)[i] - box->getAtomPosition(3)[i];}
        for (int i=0; i<3; i++) {dr13[i] = box->getAtomPosition(0)[i] - box->getAtomPosition(6)[i];}
        for (int i=0; i<3; i++) {dr23[i] = box->getAtomPosition(3)[i] - box->getAtomPosition(6)[i];}
        dist12 = sqrt(Vector::squared(dr12));
        dist13 = sqrt(Vector::squared(dr13));
        dist23 = sqrt(Vector::squared(dr23));
        dist = std::max({dist12, dist13, dist23});
        // printf("m: %d 12: %f 13: %f 23: %f max: %f \n", m, dist12, dist13, dist23, dist);
        primaryCluster.trialNotify();
        double clusterValue = primaryCluster.getValues()[0];
        printf("%f %f \n", dist, clusterValue);
      }

    }
    if (count == 3000) exit(0);
    // int iDist = int(dist);
    // if (iDist < 100) histogram[iDist]++;
    // counter++;
    // if (counter % 10000 == 0)
    // {
    //   for (int i = 0 ; i< 100; i++)
    //   {
    //     if (histogram[i] > 0) printf("%d %f\n", i, histogram[i]/counter);
    //   }





  }

  return data;
}