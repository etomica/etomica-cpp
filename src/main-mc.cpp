#include <stdio.h>

#include "potential-master.h"
#include "integrator.h"
#include "potential.h"
#include "move.h"
#include "move-volume.h"
#include "box.h"
#include "meter.h"
#include "data-sink.h"
#include "random.h"
#include "data-pump.h"
#include "util.h"
#include "action.h"

int main(int argc, char** argv) {
  int numAtoms = 500;
  double temperature = 1.0;
  double density = 1.0;
  long steps = 1000000;
  double mu = -3.5;
  double pressure = 3.7;
  bool doData = true;
  bool doHMA = true, doP = false;
  bool doGC = false, constPressure = false;
  if (doHMA) doP = false;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  PotentialLJ plj(1,1,TRUNC_SIMPLE, 3.0);
  SpeciesList speciesList;
  speciesList.add(new SpeciesSimple(1,1));
  Box box(speciesList);
  double L = pow(numAtoms/density, 1.0/3.0);
  printf("box size: %f\n", L);
  box.setBoxSize(L,L,L);
  box.setNumMolecules(0, numAtoms);
  box.initCoordinates();
  PotentialMasterCell potentialMaster(speciesList, box, false, 2);
  potentialMaster.setDoTruncationCorrection(false);
  potentialMaster.setPairPotential(0, 0, &plj);
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
  MeterPotentialEnergy meterPE(integrator);
  MCMoveVolume moveVolume(box, potentialMaster, rand, pressure, 0.1, speciesList, meterPE);
  if (constPressure) {
    integrator.addMove(&moveVolume, 1.0/numAtoms);
  }
  integrator.setTemperature(temperature);
  integrator.reset();
  PotentialCallbackHMA pcHMA(box, temperature, temperature*9.550752245164025e+00, true);
  pcHMA.setReturnAnharmonic(true, &potentialMaster);
  printf("u: %f\n", integrator.getPotentialEnergy());
  if (doData) integrator.doSteps(steps/10);
  DataPump pumpPE(meterPE, 10);
  MeterFullCompute meterFull(potentialMaster);
  PotentialCallbackPressure pcp(box, temperature, false);
  if (doP) meterFull.addCallback(&pcp);
  if (doHMA) meterFull.addCallback(&pcHMA);
  DataPump pumpFull(meterFull, 4*numAtoms);
  MeterNumAtoms meterNA(box);
  DataPump pumpNA(meterNA, 10);
  MeterDensity meterDensity(box);
  DataPump pumpDensity(meterDensity, 10);
  if (doData) {
    integrator.addListener(&pumpPE);
    if (doP || doHMA) integrator.addListener(&pumpFull);
    if (doGC) integrator.addListener(&pumpNA);
    if (constPressure) integrator.addListener(&pumpDensity);
  }

  double t1 = getTime();
  integrator.doSteps(steps);
  double t2 = getTime();
  if (doData) {
    double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
    if (!doGC) {
      statsPE[AVG_AVG] /= numAtoms;
      statsPE[AVG_ERR] /= numAtoms;
    }
    printf("u avg: %f  err: %f  cor: %f\n", statsPE[AVG_AVG], statsPE[AVG_ERR], statsPE[AVG_ACOR]);
    if (doP && !doHMA) {
      double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();
      double* statsP = statsFull[0];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);
    }
    else if (doHMA) {
      double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();
      double* statsP = statsFull[1];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);
      double* statsUHMA = statsFull[2];
      printf("uHMA avg: %f  err: %f  cor: %f\n", statsUHMA[AVG_AVG]/numAtoms, statsUHMA[AVG_ERR]/numAtoms, statsUHMA[AVG_ACOR]);
      double* statsPHMA = statsFull[3];
      double** cov = ((Average*)pumpFull.getDataSink(0))->getBlockCovariance();
      double uCor = cov[2][3]/sqrt(cov[2][2]*cov[3][3]);
      printf("pHMA avg: %f  err: %f  cor: %f   uCor: %f\n", statsPHMA[AVG_AVG], statsPHMA[AVG_ERR], statsPHMA[AVG_ACOR], uCor);

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
    if (doGC) {
      double* statsNA = ((Average*)pumpNA.getDataSink(0))->getStatistics()[0];
      printf("N avg: %f  err: %f  cor: %f\n", statsNA[AVG_AVG], statsNA[AVG_ERR], statsNA[AVG_ACOR]);
    }
    if (constPressure) {
      double* statsDensity = ((Average*)pumpDensity.getDataSink(0))->getStatistics()[0];
      printf("density avg: %f  err: %f  cor: %f\n", statsDensity[AVG_AVG], statsDensity[AVG_ERR], statsDensity[AVG_ACOR]);
    }
  }
  printf("time: %4.3f\n", t2-t1);
}
