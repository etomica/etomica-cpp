/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "data-sink.h"
#include "alloc2d.h"
#include "meter.h"
#include "vector.h"

DataWriter::DataWriter(int n, Meter* xm, const char* fn) : xMeter(xm), nData(n) {
  f = fopen(fn, "w");
}

DataWriter::~DataWriter() {
  dispose();
}

void DataWriter::dispose() {
  fclose(f);
}

void DataWriter::addData(double* x) {
  if (iLastTheta == 0 && count > 1000) return;

  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box->getMoleculeInfo(1, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  double* rfirst = box->getAtomPosition(firstAtom);
  double* rlast = box->getAtomPosition(lastAtom);
  double bvec[] = {rlast[0]-rfirst[0], rlast[1]-rfirst[1], rlast[2]-rfirst[2]};
  Vector::normalize(bvec);
  double theta = atan2(bvec[2], bvec[0]);
  int iTheta = (int)round(theta / (M_PI / 1000));
  if (iTheta == iLastTheta) return;
  iLastTheta = iTheta;
  fprintf(f, "%d %f", count, theta);
  for (int i=0; i<nData; i++)
  {
    fprintf(f, " %e", x[i]);
  }
  fprintf(f, "\n");

  count++;
}



