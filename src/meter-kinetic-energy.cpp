#include "meter.h"

MeterKineticEnergy::MeterKineticEnergy(Box& b) : Meter(1), box(b), integrator(nullptr) {
  data[0] = data[1] = 0;
}

void MeterKineticEnergy::setIntegrator(Integrator* i) {
  integrator = i;
  nData = 2;
}

double* MeterKineticEnergy::getData() {
  double imass = 1;
  double totalKE = 0;
  int numAtoms = box.getNumAtoms();
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      totalKE += imass*vi[j]*vi[j];
    }
  }
  data[0] = 0.5*totalKE;
  if (integrator) {
    double u = integrator->getPotentialEnergy();
    data[1] = data[0] + u;
  }
  return data;
}
