#include <stdio.h>

#include "potential-master.h"

PotentialCallbackEnergy::PotentialCallbackEnergy() {
  callFinished = true;
}

int PotentialCallbackEnergy::getNumData() {return 1;}

void PotentialCallbackEnergy::allComputeFinished(double uTot, double virialTot, double** f) {
  data[0] = uTot;
}

double* PotentialCallbackEnergy::getData() {
  return data;
}
