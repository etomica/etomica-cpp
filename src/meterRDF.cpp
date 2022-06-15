/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "meter.h"
#include "box.h"
#include "math.h"

MeterRDF::MeterRDF(Box& b, double r, int n) : Meter(n), box(b), rMax(r) {
  reset();
}

MeterRDF::~MeterRDF() {
  dispose();
}

void MeterRDF::dispose() {
  free(data);
  free(rData);
}

void MeterRDF::reset() {
  data = (double*)malloc(nData*sizeof(double));
  rData = (double*)malloc(nData*sizeof(double));
  setup();
}

void MeterRDF::setup() {
  for (int i=0; i<nData; i++) {
    data[i] = 0;
    rData[i] = (i+0.5)*rMax/nData;
  }
}

double* MeterRDF::getData() {
  for (int i=0; i<nData; i++) data[i] = 0;
  int na = box.getNumAtoms();
  double r2Max = rMax*rMax;
  for (int ia=0; ia<na; ia++) {
    double* ri = box.getAtomPosition(ia);
    for (int ja=ia+1; ja<na; ja++) {
      double* rj = box.getAtomPosition(ja);
      double r2 = 0;
      for (int k=0; k<3; k++) {
        double dxyz = rj[k] - ri[k];
        r2 += dxyz*dxyz;
      }
      if (r2 > r2Max) continue;
      int idx = sqrt(r2) * nData / rMax;
      data[idx]++;
    }
  }
  double density = na / box.volume();
  for (int i=0; i<nData; i++) {
    double ri = i*rMax/nData;
    double ro = (i+1)*rMax/nData;
    double Vi = 4.0/3.0*M_PI*(ro*ro*ro-ri*ri*ri);
    data[i] /= Vi;  // now we have N / V
    data[i] /= density;
  }
  return data;
}

double* MeterRDF::getRData() {
  return rData;
}

void MeterRDF::setNBins(int nb) {
  nData = nb;
  dispose();
  reset();
}

void MeterRDF::setRMax(double rm) {
  rMax = rm;
  setup();
}
