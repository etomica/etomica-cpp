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
  int numAtoms = 864;
  double temperature = 6000;
  double density = 0.15;
  long steps = 1000000;
  bool doData = true;
  bool doHMA = true;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  PotentialSSfloatTab p2(1.2446188036318708E7, 8.7932, TRUNC_SIMPLE, 6, 1000);
  PotentialSSfloatTab pRho(26068.513192447575, 8.14475, TRUNC_SIMPLE, 6, 1000);
  EmbedFsqrt embedF(6129.634374295454);
  SpeciesList speciesList;
  SpeciesSimple species(1,55.845);
  speciesList.add(&species);
  Box box(speciesList);
  double L = pow(numAtoms/density, 1.0/3.0);
  printf("box size: %f\n", L);
  box.setBoxSize(L,L,L);
  box.setNumMolecules(0, numAtoms);
  box.initCoordinates();
  box.enableVelocities();
  PotentialMasterCell potentialMaster(speciesList, box, true, 2);
  potentialMaster.setPairPotential(0, 0, &p2);
  potentialMaster.setRhoPotential(0, &pRho);
  potentialMaster.setEmbedF(0, &embedF);
  potentialMaster.init();
  IntegratorMC integrator(potentialMaster, rand);
  integrator.setTemperature(temperature);
  integrator.reset();
  MCMoveDisplacement move(box, potentialMaster, rand, 0.2);
  integrator.addMove(&move, 1);
  PotentialCallbackHMA pcHMA(box, temperature, 9.550752245164025e+00, false);
  printf("u: %f  %f\n", integrator.getPotentialEnergy()/numAtoms, integrator.getPotentialEnergy());
  if (doData) integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 1);
  MeterFullCompute meterFull(potentialMaster);
  PotentialCallbackPressure pcp(box, temperature, true);
  meterFull.addCallback(&pcp);
  if (doHMA) {
    meterFull.addCallback(&pcHMA);
  }
  DataPump pumpFull(meterFull, 4*numAtoms);
  if (doData) {
    integrator.addListener(&pumpPE);
    if (doHMA) integrator.addListener(&pumpFull);
  }

  double t1 = getTime();
  integrator.doSteps(steps);
  double t2 = getTime();
  if (doData) {
    double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
    statsPE[AVG_AVG] /= numAtoms;
    statsPE[AVG_ERR] /= numAtoms;
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
  }
  printf("time: %4.3f\n", t2-t1);
}

