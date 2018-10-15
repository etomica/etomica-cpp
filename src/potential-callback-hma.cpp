#include <stdio.h>

#include "potential-callback.h"
#include "potential-master.h"
#include "alloc2d.h"

PotentialCallbackHMA::PotentialCallbackHMA(Box& b, double T, double Ph, bool d2) : box(b), temperature(T), Pharm(Ph), phiSum(0), doD2(d2) {
  callFinished = true;
  takesForces = true;
  callPair = d2;
  data = (double*) malloc((d2 ? 6 : 4)*sizeof(double));
  int N = box.getNumAtoms();
  latticePositions = (double**)malloc2D(N, 3, sizeof(double));
  // let's pretend we can't just copy the whole thing at once
  for (int i=0; i<N; i++) {
    double* ri = box.getAtomPosition(i);
    std::copy(ri, ri+3, latticePositions[i]);
  }
}

PotentialCallbackHMA::~PotentialCallbackHMA() {
  free2D((void**)latticePositions);
  free(data);
}

int PotentialCallbackHMA::getNumData() {return doD2 ? 6 : 4;}

void PotentialCallbackHMA::pairCompute(int iAtom, int jAtom, double* drij, double u, double du, double d2u) {
  double dri[3], drj[3];
  if (phiSum == 0) {
    double *r0 = box.getAtomPosition(0);
    for (int k=0; k<3; k++) {
      dr0[k] = r0[k] - latticePositions[0][k];
    }
  }
  double* ri = box.getAtomPosition(iAtom);
  double* rj = box.getAtomPosition(jAtom);
  double r2 = 0;
  for (int k=0; k<3; k++) {
    dri[k] = ri[k] - latticePositions[iAtom][k] - dr0[k];
    drj[k] = rj[k] - latticePositions[jAtom][k] - dr0[k];
    r2 += drij[k]*drij[k];
  }
  box.nearestImage(dri);
  box.nearestImage(drj);

  double dfac = (du - d2u) / (r2*r2);
  for (int k=0; k<3; k++) {
    for (int l=0; l<3; l++) {
      double der2 = drij[k]*drij[l]*dfac;
      if (k==l) der2 -= du/r2;
      phiSum += (2*dri[k]*drj[l] - dri[k]*dri[l] - drj[k]*drj[l])*der2;
    }
  }
}

void PotentialCallbackHMA::allComputeFinished(double uTot, double virialTot, double** f) {
  const double* bs = box.getBoxSize();
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
  data[0] = uTot;
  data[1] = N*temperature/vol - virialTot/(3*vol);
  data[2] = 1.5*(N-1)*temperature + uTot + 0.5*fdotdrTot;

  double fV = (Pharm/temperature - N/vol)/(3*(N-1));
  data[3] = Pharm - virialTot/(3*vol) + fV * fdotdrTot;

  if (!doD2) return;

  double uLat = N*(-7.3212105500315605);
  double u0 = 1.5*(N-1)*temperature + uLat;
  data[4] = (uTot-u0)*(uTot-u0);
  data[5] = 1.5*(N-1) - 0.25 * (fdotdrTot + phiSum)/temperature + (data[2]-u0)*(data[2]-u0)/(temperature*temperature);

  phiSum = 0;
}

double* PotentialCallbackHMA::getData() {
  return data;
}
