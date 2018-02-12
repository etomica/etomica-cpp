#include "meter.h"

MeterPotentialEnergy::MeterPotentialEnergy(IntegratorMC& i) : Meter(1), integrator(i) { }

double* MeterPotentialEnergy::getData() {
  data[0] = integrator.getPotentialEnergy();
  return data;
}
