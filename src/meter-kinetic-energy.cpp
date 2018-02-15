#include "meter.h"

MeterKineticEnergy::MeterKineticEnergy(Box& b) : Meter(1), box(b) { }

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
  return data;
}
