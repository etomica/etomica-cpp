#include "meter.h"
#include "box.h"

MeterNumAtoms::MeterNumAtoms(Box& b) : Meter(1), box(b) { }

double* MeterNumAtoms::getData() {
  data[0] = box.getNumAtoms();
  return data;
}
