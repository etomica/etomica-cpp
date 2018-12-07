/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>

#include "potential-callback.h"
#include "box.h"

PotentialCallbackVirialTensor::PotentialCallbackVirialTensor(Box& b) : box(b) {
  callFinished = true;
  takesForces = true;
  takesVirialTensor = true;
}

int PotentialCallbackVirialTensor::getNumData() {return 6;}

void PotentialCallbackVirialTensor::allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor) {
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  
  for (int i=0; i<6; i++) {
    data[i] = virialTensor[i]/vol;
  }
}

double* PotentialCallbackVirialTensor::getData() {
  return data;
}
