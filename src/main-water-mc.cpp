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
#include "action.h"

int main(int argc, char** argv) {
  int numMolecules = 46;
  double K = 0.8314459861448581;
  double temperature = 120*K;
  long steps = 10000000;
  bool doData = true;
  bool doHMA = true;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  double A = 600E3, C = 610;
  double s6 = A/C;
  double sigma = pow(s6, 1.0/6.0);
  double epsilon = (C/s6)*1000/4 * 1000.*1e20*1e-24*4.184;
  printf("sigma: %f\n", sigma);
  printf("epsilon: %f\n", epsilon);
  PotentialLJ pOO(epsilon,sigma,TRUNC_SIMPLE, 11);
  double alpha = 0.26111648393354675;
  double kCut = 1.5;
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
  box.enableVelocities();
  PotentialMasterCell potentialMaster(speciesList, box, false, 3);
  potentialMaster.setPairPotential(oType, oType, &pOO);
  potentialMaster.setPairPotential(hType, hType, &pHH);
  potentialMaster.setPairPotential(hType, mType, &pHM);
  potentialMaster.setPairPotential(mType, mType, &pMM);
  potentialMaster.setCharge(hType, qH);
  potentialMaster.setCharge(mType, -2*qH);
  potentialMaster.setEwald(kCut, alpha);
  potentialMaster.setDoTruncationCorrection(false);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);
  PotentialCallbackEnergy pce;
  /*MeterFullCompute meterFull0(potentialMaster);
  meterFull0.setDoCompute(true);
  meterFull0.addCallback(&pce);
  for (int i=0; i<1000; i++) {
    meterFull0.getData();
  }
  double* data = meterFull0.getData();
  printf("u %f %f\n", data[0], data[0]/numMolecules);*/
  IntegratorMC integrator(potentialMaster, rand);
  //IntegratorNHC integrator(speciesList.getAtomInfo(), potentialMaster, rand, box, 3, 0.1);
  //integrator.setTimeStep(0.005);
  integrator.setTemperature(temperature);
  integrator.reset();
  MCMoveMoleculeDisplacement move(box, potentialMaster, rand, 0.2);
  integrator.addMove(&move, 1);
  MCMoveMoleculeRotate moveRotate(speciesList, box, potentialMaster, rand);
  integrator.addMove(&moveRotate, 1);
  PotentialCallbackMoleculeHMA pcHMA(box, speciesList, temperature);
  double u0 = integrator.getPotentialEnergy();
  printf("u: %f\n", u0/numMolecules);
  if (doData) integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 1);
  MeterFullCompute meterFull(potentialMaster);
  meterFull.setDoCompute(true);
  if (doHMA) {
    meterFull.addCallback(&pcHMA);
  }
  DataPump pumpFull(meterFull, 4*numMolecules);
  if (doData) {
    integrator.addListener(&pumpPE);
    if (doHMA) integrator.addListener(&pumpFull);
  }

  double t1 = getTime();
  integrator.doSteps(steps);
  double t2 = getTime();
  if (doData) {
    double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
    statsPE[AVG_AVG] /= numMolecules;
    statsPE[AVG_ERR] /= numMolecules;
    printf("u avg: %f  err: %f  cor: %f\n", statsPE[AVG_AVG], statsPE[AVG_ERR], statsPE[AVG_ACOR]);
    if (doHMA) {
      double uah = (statsPE[AVG_AVG] * numMolecules - u0 - (1.5*(numMolecules-1) + 1.5*numMolecules)*temperature) / numMolecules;
     
      printf("uah avg: %f  err: %f  cor: %f\n", uah, statsPE[AVG_ERR], statsPE[AVG_ACOR]);
      double* statsHMA = ((Average*)pumpFull.getDataSink(0))->getStatistics()[0];
      statsHMA[AVG_AVG] /= numMolecules;
      statsHMA[AVG_ERR] /= numMolecules;
      printf("uHMA avg: %f  err: %f  cor: %f\n", statsHMA[AVG_AVG], statsHMA[AVG_ERR], statsHMA[AVG_ACOR]);
      uah = (statsHMA[AVG_AVG] * numMolecules - u0 - (1.5*(numMolecules-1) + 1.5*numMolecules)*temperature) / numMolecules;
      printf("uahHMA avg: %f  err: %f  cor: %f\n", uah, statsHMA[AVG_ERR], statsHMA[AVG_ACOR]);
    }
  }
  printf("time: %4.3f\n", t2-t1);
}

