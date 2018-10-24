/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "meter.h"
#include "box.h"

MeterDensity::MeterDensity(Box& b) : Meter(1), box(b) { }

double* MeterDensity::getData() {
  const double *bs = box.getBoxSize();
  data[0] = box.getNumAtoms()/(bs[0]*bs[1]*bs[2]);
  return data;
}
