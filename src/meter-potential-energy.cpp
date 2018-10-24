/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "meter.h"

MeterPotentialEnergy::MeterPotentialEnergy(Integrator& i) : Meter(1), integrator(i) { }

double* MeterPotentialEnergy::getData() {
  data[0] = integrator.getPotentialEnergy();
  return data;
}
