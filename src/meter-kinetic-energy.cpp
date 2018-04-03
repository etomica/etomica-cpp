#include "meter.h"

MeterKineticEnergy::MeterKineticEnergy(Box& b) : Meter(1), box(b), integrator(nullptr) {
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
