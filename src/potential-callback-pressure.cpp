#include <stdio.h>

#include "potential-master.h"

PotentialCallbackPressure::PotentialCallbackPressure(Box& b, double T) : box(b), temperature(T) {
  callFinished = true;
}

int PotentialCallbackPressure::getNumData() {return 1;}

void PotentialCallbackPressure::allComputeFinished(double uTot, double virialTot, double** f) {
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  int numAtoms = box.getNumAtoms();
  
  data[0] = (numAtoms*temperature - virialTot/3)/vol;
}

double* PotentialCallbackPressure::getData() {
  return data;
}
