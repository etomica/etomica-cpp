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
  double temperature = 1;
  double density = 1.0;
  long steps = 20000;
  bool doData = true;
  bool doHMA = false, doHMA2 = false;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  PotentialLJ plj(1,1,TRUNC_FORCE_SHIFT, 3.0);
  SpeciesList speciesList;
  speciesList.add(new SpeciesSimple(1,1));
  Box box(speciesList);
  double L = pow(numAtoms/density, 1.0/3.0);
  printf("box size: %f\n", L);
  box.setBoxSize(L,L,L);
  box.setNumMolecules(0, numAtoms);
  box.initCoordinates();
  box.enableVelocities();
  /*PotentialMasterCell potentialMaster(plj, box, 3.0, 2);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);*/
  PotentialMasterList potentialMaster(speciesList, box, false, 2, 3.5);
  potentialMaster.setDoTruncationCorrection(false);
  potentialMaster.setPairPotential(0, 0, &plj);
  potentialMaster.init();
  //PotentialMaster potentialMaster(plj, box);
  //IntegratorNVE integrator(speciesList.getAtomInfo(), potentialMaster, rand, box);
  IntegratorNHC integrator(speciesList.getAtomInfo(), potentialMaster, rand, box, 3, 0.1);
  integrator.setTimeStep(0.005);
  integrator.setTemperature(temperature);
  integrator.setNbrCheckInterval(20);
  integrator.randomizeVelocities(true);
  integrator.reset();
  PotentialCallbackHMA pcHMA(box, temperature, 9, doHMA2);
  double uLat = integrator.getPotentialEnergy();
  printf("u: %f\n", uLat);
  if (doData) integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 1);
  MeterFullCompute meterFull(potentialMaster);
  meterFull.setDoCompute(false);
  PotentialCallbackPressure pcp(box, temperature, true);
  if (doHMA) {
    meterFull.addCallback(&pcHMA);
  }
  else {
    meterFull.addCallback(&pcp);
  }
  DataPump pumpFull(meterFull, 10);
  MeterKineticEnergy meterKE;
  meterKE.setIntegrator(&integrator);
  DataPump pumpKE(meterKE, 10);
  if (doData) {
    integrator.addListener(&pumpPE);
    integrator.addListener(&pumpFull);
    if (doHMA) {
      integrator.addPotentialCallback(&pcHMA, 10);
    }
    else {
      integrator.addPotentialCallback(&pcp);
    }

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
    if (!doHMA) {
      double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();
      double* statsP = statsFull[0];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);
    }
    else {
      double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();
      double* statsP = statsFull[1];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);
      if (doHMA2) {
        double* statsU = statsFull[0];
        double* statsU2 = statsFull[4];
        double uvar = statsU2[AVG_AVG] - statsU[AVG_AVG]*statsU[AVG_AVG];
        double Cv = uvar/(temperature*temperature);
        double** cor = ((Average*)pumpFull.getDataSink(0))->getBlockCorrelation();
        double x = statsU[AVG_AVG];
        double ex = statsU[AVG_ERR];
        double eCv = sqrt(statsU2[AVG_ERR]*statsU2[AVG_ERR] + 4*x*x*ex*ex
                        - 4 * x*ex*statsU2[AVG_ERR]*cor[0][4])/(temperature*temperature);
        printf("Cv avg: %f  err: %f  cor: %f\n", Cv/numAtoms, eCv/numAtoms, statsU2[AVG_ACOR]);
      }
      double* statsUHMA = statsFull[2];
      printf("uHMA avg: %17.10e  err: %10.4e  cor: %f\n", statsUHMA[AVG_AVG]/numAtoms, statsUHMA[AVG_ERR]/numAtoms, statsUHMA[AVG_ACOR]);
      double* statsPHMA = statsFull[3];
      printf("pHMA avg: %17.10e  err: %10.4e  cor: %f\n", statsPHMA[AVG_AVG], statsPHMA[AVG_ERR], statsPHMA[AVG_ACOR]);
      if (doHMA2) {
        double* statsCvHMA = statsFull[5];
        printf("CvHMAraw: %20.15e  err: %10.4e  cor: %f\n", statsCvHMA[AVG_AVG]/numAtoms, statsCvHMA[AVG_ERR]/numAtoms, statsCvHMA[AVG_ACOR]);
        double x = statsUHMA[AVG_AVG]/temperature;
        double Cv = statsCvHMA[AVG_AVG] - x*x;
        double** cor = ((Average*)pumpFull.getDataSink(0))->getBlockCorrelation();
        double ex = statsUHMA[AVG_ERR]/temperature;
        double eCv = sqrt(statsCvHMA[AVG_ERR]*statsCvHMA[AVG_ERR] + 4*x*x*ex*ex
                        - 4 * x*ex*statsCvHMA[AVG_ERR]*cor[2][5]);
        printf("CvHMA avg: %17.10e  err: %10.4e  cor: %f\n", Cv/numAtoms, eCv/numAtoms, statsCvHMA[AVG_ACOR]);
      }
    }
  }
  printf("time: %4.3f\n", t2-t1);
}

