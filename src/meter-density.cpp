#include "meter.h"

MeterDensity::MeterDensity(Box& b) : Meter(1), box(b) { }

double* MeterDensity::getData() {
  const double *bs = box.getBoxSize();
  data[0] = box.getNumAtoms()/(bs[0]*bs[1]*bs[2]);
  return data;
}
