/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>

#include "potential-master.h"

PotentialCallbackEnergy::PotentialCallbackEnergy() {
  callFinished = true;
}

int PotentialCallbackEnergy::getNumData() {return 1;}

void PotentialCallbackEnergy::allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor) {
  data[0] = uTot;
}

double* PotentialCallbackEnergy::getData() {
  return data;
}
