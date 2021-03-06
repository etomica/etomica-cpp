/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>

#include "potential-master.h"

PotentialCallbackPressure::PotentialCallbackPressure(Box& b, double T, bool tf) : box(b), temperature(T) {
  callFinished = true;
  takesForces = tf;
}

int PotentialCallbackPressure::getNumData() {return 1;}

void PotentialCallbackPressure::allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor) {
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  int numMolecules = box.getTotalNumMolecules();
  
  data[0] = (numMolecules*temperature - virialTot/3)/vol;
}

double* PotentialCallbackPressure::getData() {
  return data;
}
