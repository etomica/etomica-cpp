/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "box.h"
#include "ewald.h"
#include "species.h"
#include "alloc2d.h"
#include "util.h"

EwaldBase::EwaldBase(const SpeciesList &sl, Box& b) : speciesList(sl), box(b), rigidMolecules(true) {
  int numAtomTypes = sl.getNumAtomTypes();
  charges = new double[numAtomTypes];
  for (int i=0; i<numAtomTypes; i++) charges[i] = 0;
  B6 = (double**)malloc2D(numAtomTypes, numAtomTypes, sizeof(double));
  b6 = (double**)malloc2D(numAtomTypes, 6, sizeof(double));
  for (int i=0; i<numAtomTypes; i++) {
    for (int j=0; j<numAtomTypes; j++) B6[i][j] = 0;
    for (int j=0; j<6; j++) b6[i][j] = 0;
  }
}

EwaldBase::~EwaldBase() {
  delete[] charges;
  free2D((void**)B6);
  free2D((void**)b6);
}

void EwaldBase::setCharge(int iType, double q) {
  charges[iType] = q;
}

void EwaldBase::setR6Coeff(int iType, double sigma, double epsilon) {
  int numAtomTypes = speciesList.getNumAtomTypes();
  double sigmak = 1;
  for (int jType=0; jType<numAtomTypes; jType++) {
    B6[iType][jType] = 0;
    for (int k=0; k<6; k++) {
      int ck = factorial(6)/(factorial(6-k)*factorial(k));
      b6[iType][k] = 0.25*sigmak*sqrt(ck*epsilon);
      B6[iType][jType] += b6[iType][k]*b6[jType][6-k];
      sigmak *= sigma;
    }
    B6[jType][iType] = B6[iType][jType];
  }
}

void EwaldBase::setRigidMolecules(bool rm) {
  rigidMolecules = rm;
}
