/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>

#include "potential-master.h"
#include "integrator.h"
#include "potential.h"
#include "move.h"
#include "box.h"
#include "meter.h"
#include "util.h"
#include "lattice-dynamics-full.h"
#include "action.h"
#include "ewald.h"
#include "alloc2d.h"

void usage() {
  fprintf(stderr, "usage: ld-fi [--ss|--ss6] [--shift|--fshift] [--atomicExp] ncells density rc nwv\n");
}

int main(int argc, char** argv) {
  if (argc < 5 || argc > 8) {
    usage();
    exit(1);
  }
  int arg0 = 1;
  bool ss = false, ss6 = false;
  if (!strcmp(argv[1], "--ss")) {
    arg0++;
    ss = true;
  }
  else if (!strcmp(argv[1], "--ss6")) {
    arg0++;
    ss6 = true;
  }
  int trunc = TRUNC_SIMPLE;
  if (!strcmp(argv[arg0], "--shift")) {
    trunc = TRUNC_SHIFT;
    printf(" shifted\n");
    arg0++;
  }
  else if (!strcmp(argv[arg0], "--fshift")) {
    trunc = TRUNC_FORCE_SHIFT;
    printf(" force-shifted\n");
    arg0++;
  }

  bool atomicExp = false;
  if (!strcmp(argv[arg0], "--atomicExp")) {
    arg0++;
    atomicExp = true;
    printf(" doing atomic exponential: exp(-i k.r_ij) \n");
  }

  int nCells = atoi(argv[arg0]);
  char* c2 = argv[arg0+1];
  char* c1 = c2;
  double density = strtod(c1, &c2);
  if (c1==c2) {
    usage();
    exit(1);
  }
  c2 = argv[arg0+2];
  c1 = c2;
  double rc = strtod(c1, &c2);
  if (c1==c2) {
    usage();
    exit(1);
  }
  int nwv = atoi(argv[arg0+3]);
  Potential* p = nullptr;
  if (ss) p = new PotentialSS(4.0, 12, trunc, rc);
  else if (ss6) p = new PotentialSS(4.0, 6, trunc, rc);
  else p = new PotentialLJ(1.0, 1.0, trunc, rc);
  int nBasis = 4;
  int numAtoms = nBasis*pow(nCells, 3);
  SpeciesList speciesList;
  speciesList.add(new SpeciesSimple(1,1));
  Box box(speciesList);
  double cellSize = pow(nBasis/density, 1.0/3.0);
  double L = cellSize*nCells;
  box.setBoxSize(L,L,L);
  box.setNumMolecules(0, numAtoms);
  box.initCoordinates();
  PotentialMasterCell potentialMaster(speciesList, box, false, 2);
  potentialMaster.setDoTruncationCorrection(false);
  potentialMaster.setPairPotential(0, 0, p);
  potentialMaster.init();
  //int* numCells = potentialMaster.getNumCells();
  //printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);
  PotentialCallbackEnergy pce;
  MeterFullCompute meterFull(potentialMaster);
  meterFull.addCallback(&pce);
  double* data = meterFull.getData();
  double lastU = data[0]/numAtoms;

  double wvMax = 2.0*M_PI / cellSize;
  long nwv_tot = pow(nwv, 3);
  double** waveVectors = (double**)malloc2D(nwv_tot, 3, sizeof(double));
  int l = 0;
  for (int i=0; i<nwv; i++) {
    for (int j=0; j<nwv; j++) {
      for (int k=0; k<nwv; k++) {
        waveVectors[l][0] = i*wvMax/nwv;
        waveVectors[l][1] = j*wvMax/nwv;
        waveVectors[l][2] = k*wvMax/nwv;
        l++;
      }
    }
  }

  LatticeDynamicsFull ld(potentialMaster, atomicExp);
  ld.setNumCells(nCells, nCells, nCells);
  ld.setupForWV(nwv_tot, waveVectors);

  vector<PotentialCallback*> pcs;
  pcs.push_back(&ld);
  potentialMaster.computeAll(pcs);



  //double t1 = getTime();
  //ld.compute();
  //double t2 = getTime();

  if  (ld.getUnstable()) {
    printf("System is unstable\n");
    exit(1);
  }

  long N_eff = nwv_tot*nBasis;
  double A = ld.getLogSum()/2.0;
  printf(" N_eff  %ld \n", N_eff);
  printf(" nc  %d \n", nCells);
  printf(" a  %15.10f \n", cellSize);
  printf(" u0 %22.15f\n", lastU);
  printf(" A/N_eff  %22.15f\n", A/N_eff);
  printf(" A/N_eff-ln(N_eff)/N_eff  %22.15f\n", A/N_eff - log(N_eff)/N_eff);

//  printf("time %4.3f\n", t2-t1);
}
