#include <stdio.h>

#include "potential-master.h"
#include "integrator.h"
#include "potential.h"
#include "move.h"
#include "box.h"
#include "meter.h"
#include "data-sink.h"
#include "random.h"
#include "util.h"

int main(int argc, char** argv) {
  int numAtoms = 864;
  double temperature = 8000;
  double density = 0.15;
  long steps = 1000;
  bool doData = true;
  bool doHMA = false;

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
  PotentialMasterList potentialMaster(speciesList, box, true, 2, 7.2);
  potentialMaster.setPairPotential(0, 0, &p2);
  potentialMaster.setRhoPotential(0, &pRho);
  potentialMaster.setEmbedF(0, &embedF);
  potentialMaster.init();
  //PotentialMaster potentialMaster(plj, box);
  IntegratorMD integrator(speciesList.getAtomInfo(), potentialMaster, rand, box);
  integrator.setTimeStep(0.001);
  integrator.setTemperature(temperature);
  integrator.setNbrCheckInterval(1);
  integrator.reset();
  PotentialCallbackHMA pcHMA(box, temperature, 9.550752245164025e+00);
  printf("u: %f\n", integrator.getPotentialEnergy()/numAtoms);
  if (doData) integrator.doSteps(steps/10);
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
  DataPump pumpKE(meterKE, 10);
  if (doData) {
    integrator.addListener(&pumpPE);
    if (doHMA) integrator.addListener(&pumpFull);
    integrator.addListener(&pumpKE);
  }

  double t1 = getTime();
  integrator.doSteps(steps);
  double t2 = getTime();
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
  printf("time: %4.3f\n", t2-t1);
}

