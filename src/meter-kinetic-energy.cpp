/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "meter.h"

MeterKineticEnergy::MeterKineticEnergy() : Meter(1), integrator(nullptr) {
  data[0] = data[1] = 0;
}

void MeterKineticEnergy::setIntegrator(IntegratorMD* i) {
  integrator = i;
  nData = 2;
}

double* MeterKineticEnergy::getData() {
  data[0] = integrator->getKineticEnergy();
  if (integrator) {
    double u = integrator->getPotentialEnergy();
    data[1] = data[0] + u;
  }
  return data;
}
