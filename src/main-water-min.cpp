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
#include "minimize.h"
#include "action.h"
#include "ewald.h"

int main(int argc, char** argv) {
  int numMolecules = 46;

  double A = 600E3, C = 610;
  double s6 = A/C;
  double sigma = pow(s6, 1.0/6.0);
  double epsilon = (C/s6)*1000/4 * 1000.*1e20*1e-24*4.184;
  PotentialLJ pOO(epsilon,sigma,TRUNC_SIMPLE, 20);
  double alpha = 0.3;
  double kCut = 3.5;
  double qH = 193.82504408037946;
  PotentialEwaldBare pHH(alpha, qH*qH, 11);
  PotentialEwaldBare pHM(alpha, -2*qH*qH, 11);
  PotentialEwaldBare pMM(alpha, 4*qH*qH, 11);
  SpeciesList speciesList;
  SpeciesFile species("water.species");
  speciesList.add(&species);
  int oType = species.getTypeForSymbol("O");
  int mType = species.getTypeForSymbol("M");
  int hType = species.getTypeForSymbol("H");
  Box box(speciesList);
  double L = 12.03;
  printf("box size: %f\n", L);
  box.setBoxSize(L,L,L);
  box.setNumMolecules(0, numMolecules);
  ConfigurationFile config(box, "water.pos");
  config.go();
  PotentialMasterCell potentialMaster(speciesList, box, false, 3);
  potentialMaster.setPairPotential(oType, oType, &pOO);
  potentialMaster.setPairPotential(hType, hType, &pHH);
  potentialMaster.setPairPotential(hType, mType, &pHM);
  potentialMaster.setPairPotential(mType, mType, &pMM);
  EwaldFourier ewald(speciesList, box);
  ewald.setCharge(hType, qH);
  ewald.setCharge(mType, -2*qH);
  ewald.setParameters(kCut, alpha);
  potentialMaster.setEwald(&ewald);
  potentialMaster.setDoTruncationCorrection(true);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);
  PotentialCallbackEnergy pce;
  MeterFullCompute meterFull(potentialMaster);
  meterFull.addCallback(&pce);
  double* data = meterFull.getData();
  double lastU = data[0]/numMolecules;
  printf("u0: %f\n", lastU);

  double t1 = getTime();
  Minimize min(potentialMaster);
  for (int i=0; i<3; i++) {
    min.doStep();
    double lastDR = min.getLastDR();
    data = meterFull.getData();
    double unow = data[0]/numMolecules;
    printf("%d % 10.4e % 10.4e\n", i, lastDR, unow-lastU);
    lastU = unow;
  }
  double t2 = getTime();
  printf("time: %4.3f\n", t2-t1);
}

