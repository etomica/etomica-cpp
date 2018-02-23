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
  int numAtoms = 32000;
  double temperature = 1.44;
  double density = 0.8442;
  long steps = 1000;
  bool doData = true;
  bool doHMA = false;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  PotentialLJ plj(TRUNC_SIMPLE, 2.5);
  SpeciesList speciesList;
  speciesList.add(new SpeciesSimple(1,1));
  Box box(speciesList);
  double L = pow(numAtoms/density, 1.0/3.0);
  printf("box size: %f\n", L);
  box.setBoxSize(L,L,L);
  box.boxSizeUpdated();
  box.setNumMolecules(0, numAtoms);
  box.initCoordinates();
  box.enableVelocities();
  /*PotentialMasterCell potentialMaster(plj, box, 3.0, 2);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);*/
  PotentialMasterList potentialMaster(speciesList, box, 2, 2.8);
  potentialMaster.setPairPotential(0, 0, &plj);
  potentialMaster.init();
  //PotentialMaster potentialMaster(plj, box);
  IntegratorMD integrator(potentialMaster, rand, box);
  integrator.setTimeStep(0.005);
  integrator.setTemperature(temperature);
  integrator.setNbrCheckInterval(20);
  integrator.reset();
  PotentialCallbackHMA pcHMA(box, temperature, 9.550752245164025e+00);
  printf("u: %f\n", integrator.getPotentialEnergy());
  //if (doData) integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 1);
  MeterFullCompute meterFull(potentialMaster);
  meterFull.setDoCompute(false);
  PotentialCallbackPressure pcp(box, temperature);
  meterFull.addCallback(&pcp);
  integrator.addPotentialCallback(&pcp);
  if (doHMA) {
    meterFull.addCallback(&pcHMA);
    integrator.addPotentialCallback(&pcHMA);
  }
  DataPump pumpFull(meterFull, 1);
  MeterKineticEnergy meterKE(box);
  meterKE.setIntegrator(&integrator);
  double* dataKE0 = meterKE.getData();
  printf("T0: %f\n", dataKE0[0]/(1.5*(numAtoms-1)));
  DataPump pumpKE(meterKE, 10);
  if (doData) {
    integrator.addListener(&pumpPE);
    if (doHMA) integrator.addListener(&pumpFull);
    integrator.addListener(&pumpKE);
  }

  integrator.doSteps(steps);
  if (doData) {
    double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
    statsPE[AVG_AVG] /= numAtoms;
    statsPE[AVG_ERR] /= numAtoms;
    printf("u avg: %f  err: %f  cor: %f\n", statsPE[AVG_AVG], statsPE[AVG_ERR], statsPE[AVG_ACOR]);
    double** statsKEE = ((Average*)pumpKE.getDataSink(0))->getStatistics();
    double* statsKE = statsKEE[0];
    statsKE[AVG_AVG] /= 1.5*(numAtoms-1);
    statsKE[AVG_ERR] /= 1.5*(numAtoms-1);
    printf("T avg: %f  err: %f  cor: %f\n", statsKE[AVG_AVG], statsKE[AVG_ERR], statsKE[AVG_ACOR]);
    double* statsE = statsKEE[1];
    statsE[AVG_AVG] /= numAtoms;
    statsE[AVG_ERR] /= numAtoms;
    printf("E avg: %f  err: %f  cor: %f\n", statsE[AVG_AVG], statsE[AVG_ERR], statsE[AVG_ACOR]);
    if (doHMA) {
      double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();
      double* statsP = statsFull[0];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);
      double* statsUHMA = statsFull[1];
      printf("uHMA avg: %f  err: %f  cor: %f\n", statsUHMA[AVG_AVG]/numAtoms, statsUHMA[AVG_ERR]/numAtoms, statsUHMA[AVG_ACOR]);
      double* statsPHMA = statsFull[2];
      printf("pHMA avg: %f  err: %f  cor: %f\n", statsPHMA[AVG_AVG], statsPHMA[AVG_ERR], statsPHMA[AVG_ACOR]);
    }
  }
}

