/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "meter.h"
#include "box.h"

MeterNumAtoms::MeterNumAtoms(Box& b) : Meter(1), box(b) { }

double* MeterNumAtoms::getData() {
  data[0] = box.getNumAtoms();
  return data;
}
