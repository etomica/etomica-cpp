#include <stdio.h>

#include "potential-master.h"
#include "alloc2d.h"

PotentialCallbackHMA::PotentialCallbackHMA(Box& b, double T, double Ph) : box(b), temperature(T), Pharm(Ph) {
  callFinished = true;
  takesForces = true;
  int N = box.getNumAtoms();
  latticePositions = (double**)malloc2D(N, 3, sizeof(double));
  // let's pretend we can't just copy the whole thing at once
  for (int i=0; i<N; i++) {
    double* ri = box.getAtomPosition(i);
    std::copy(ri, ri+3, latticePositions[i]);
  }
}

int PotentialCallbackHMA::getNumData() {return 2;}

void PotentialCallbackHMA::allComputeFinished(double uTot, double virialTot, double** f) {
  double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  int N = box.getNumAtoms();
  
  double dr0[3], dr[3];
  double *r0 = box.getAtomPosition(0);
  for (int j=0; j<3; j++) {
    dr0[j] = r0[j] - latticePositions[0][j];
  }
  double fdotdrTot = 0;
  for (int i=0; i<N; i++) {
    double* ri = box.getAtomPosition(i);
    for (int j=0; j<3; j++) {
      dr[j] = ri[j] - latticePositions[i][j] - dr0[j];
    }
    box.nearestImage(dr);
    for (int j=0; j<3; j++) {
      fdotdrTot += f[i][j]*dr[j];
    }
  }
  data[0] = 1.5*(N-1)*temperature + uTot + 0.5*fdotdrTot;

  double fV = (Pharm/temperature - N/vol)/(3*(N-1));
  data[1] = Pharm - virialTot/(3*vol) + fV * fdotdrTot;
}

double* PotentialCallbackHMA::getData() {
  return data;
}
