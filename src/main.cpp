#include <stdio.h>

#include "potential-master.h"
#include "integrator.h"
#include "potential.h"
#include "move.h"
#include "box.h"
#include "meter.h"
#include "data-sink.h"
#include "random.h"

int main(int argc, char** argv) {
  int numAtoms = 500;
  double temperature = 1.0;
  double density = 1.0;
  long steps = 1000000;
  double mu = -3.5;
  bool doData = false;
  bool doHMA = false;
  bool doGC = false;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  PotentialLJ plj(TRUNC_SIMPLE, 3.0);
  SpeciesList speciesList;
  speciesList.add(new SpeciesSimple(1,1));
  Box box(speciesList);
  double L = pow(numAtoms/density, 1.0/3.0);
  printf("box size: %f\n", L);
  box.setBoxSize(L,L,L);
  box.boxSizeUpdated();
  box.setNumMolecules(0, numAtoms);
  box.initCoordinates();
  PotentialMasterCell potentialMaster(plj, box, 3.0, 2);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);
  //PotentialMaster potentialMaster(plj, box);
  IntegratorMC integrator(potentialMaster, rand);
  MCMoveDisplacement move(box, potentialMaster, rand, 0.2);
  integrator.addMove(&move, 1);
  MCMoveInsertDelete moveID(box, potentialMaster, rand, mu, 0);
  if (doGC) {
    integrator.addMove(&moveID, 1);
  }
  integrator.setTemperature(temperature);
  integrator.reset();
  //PotentialCallbackHMA pcHMA(box, temperature, 9.550752245164025e+00);
  printf("u: %f\n", integrator.getPotentialEnergy());
  if (doData) integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 10);
  MeterFullCompute meterFull(potentialMaster);
  PotentialCallbackPressure pcp(box, temperature);
  meterFull.addCallback(&pcp);
  //if (doHMA) meterFull.addCallback(&pcHMA);
  DataPump pumpFull(meterFull, 4*numAtoms);
  MeterNumAtoms meterNA(box);
  DataPump pumpNA(meterNA, 10);
  if (doData) {
    integrator.addListener(&pumpPE);
    if (doHMA) integrator.addListener(&pumpFull);
    if (doGC) integrator.addListener(&pumpNA);
  }

  integrator.doSteps(steps);
  if (doData) {
    double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
    if (!doGC) {
      statsPE[AVG_AVG] /= numAtoms;
      statsPE[AVG_ERR] /= numAtoms;
    }
    printf("u avg: %f  err: %f  cor: %f\n", statsPE[AVG_AVG], statsPE[AVG_ERR], statsPE[AVG_ACOR]);
    if (doHMA) {
      double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();
      double* statsP = statsFull[0];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);
      double* statsUHMA = statsFull[1];
      printf("uHMA avg: %f  err: %f  cor: %f\n", statsUHMA[AVG_AVG]/numAtoms, statsUHMA[AVG_ERR]/numAtoms, statsUHMA[AVG_ACOR]);
      double* statsPHMA = statsFull[2];
      printf("pHMA avg: %f  err: %f  cor: %f\n", statsPHMA[AVG_AVG], statsPHMA[AVG_ERR], statsPHMA[AVG_ACOR]);
    }
    if (doGC) {
      double* statsNA = ((Average*)pumpNA.getDataSink(0))->getStatistics()[0];
      printf("N avg: %f  err: %f  cor: %f\n", statsNA[AVG_AVG], statsNA[AVG_ERR], statsNA[AVG_ACOR]);
    }
  }
}
