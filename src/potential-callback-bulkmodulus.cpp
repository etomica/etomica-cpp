/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>

#include "potential-callback.h"
#include "box.h"

/**
 * This class computes contributions to the bulk modulus.  The formula is
 * B = Var[y0]/kT +  y1
 *
 * This class returns both y0 and y1 as separate peices of data.
 */
PotentialCallbackBulkModulus::PotentialCallbackBulkModulus(Box& b, double T) : box(b), temperature(T) {
  takesForces = true;
  callFinished = true;
  callPair = true;
}

int PotentialCallbackBulkModulus::getNumData() {return 2;}

void PotentialCallbackBulkModulus::reset() {data[0] = data[1] = 0;}

void PotentialCallbackBulkModulus::pairCompute(int iAtom, int jAtom, double* drij, double u, double du, double d2u) {
  data[1] -= d2u;
}

void PotentialCallbackBulkModulus::allComputeFinished(double uTot, double virialTot, double** f) {
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  int numMolecules = box.getTotalNumMolecules();
  
  // the variance of data[0] must be subtracted from data[1]
  data[0] = -virialTot/(3*vol);
  data[1] = -(numMolecules*temperature + data[1]/9 + 2*virialTot/9)/vol;
}

double* PotentialCallbackBulkModulus::getData() {
  return data;
}
