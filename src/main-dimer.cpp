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
#include "data-sink.h"
#include "data-pump.h"
#include "random.h"
#include "util.h"

int main(int argc, char** argv) {
  int numMolecules = 500;
  double temperature = 2.0;
  double density = 0.5;
  long steps = 1000000;
  double mu = -3.5;
  bool doData = true;
  bool doGC = false;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  PotentialLJ plj(1,1,TRUNC_SIMPLE, 3.0);
  SpeciesList speciesList;
  SpeciesSimple dimer(2,1);
  dimer.setAtomPosition(0, -0.25, 0, 0);
  dimer.setAtomPosition(1, +0.25, 0, 0);
  speciesList.add(&dimer);
  Box box(speciesList);
  double L = pow(numMolecules/density, 1.0/3.0);
  printf("box size: %f\n", L);
  box.setBoxSize(L,L,L);
  box.setNumMolecules(0, numMolecules);
  box.initCoordinates();
  PotentialMasterCell potentialMaster(speciesList, box, false, 2);
  potentialMaster.setPairPotential(0, 0, &plj);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);
  //PotentialMaster potentialMaster(plj, box);
  IntegratorMC integrator(potentialMaster, rand);
  MCMoveMoleculeDisplacement move(box, potentialMaster, rand, 0.2);
  integrator.addMove(&move, 1);
  MCMoveMoleculeRotate moveRotate(speciesList, box, potentialMaster, rand);
  integrator.addMove(&moveRotate, 1);
  MCMoveInsertDelete moveID(box, potentialMaster, rand, mu, 0);
  if (doGC) {
    integrator.addMove(&moveID, 1);
  }
  integrator.setTemperature(temperature);
  integrator.reset();
  printf("u: %f\n", integrator.getPotentialEnergy());
  if (doData) integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 10);
  MeterFullCompute meterFull(potentialMaster);
  PotentialCallbackPressure pcp(box, temperature, false);
  meterFull.addCallback(&pcp);
  DataPump pumpFull(meterFull, 4*numMolecules);
  MeterNumAtoms meterNA(box);
  DataPump pumpNA(meterNA, 10);
  if (doData) {
    integrator.addListener(&pumpPE);
  }

  double t1 = getTime();
  integrator.doSteps(steps);
  double t2 = getTime();
  if (doData) {
    double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
    statsPE[AVG_AVG] /= numMolecules;
    statsPE[AVG_ERR] /= numMolecules;
    printf("u avg: %f  err: %f  cor: %f\n", statsPE[AVG_AVG], statsPE[AVG_ERR], statsPE[AVG_ACOR]);
  }
  printf("time: %4.3f\n", t2-t1);
}
