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
  double temperature = 10;
  long steps = 2000;
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
  PotentialLJ pOO(epsilon,sigma,TRUNC_SIMPLE, 12);
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
  /*PotentialMasterCell potentialMaster(plj, box, 3.0, 2);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);*/
  PotentialMasterList potentialMaster(speciesList, box, false, 2, 12);
  printf("O %d  H %d  M %d\n", oType, hType, mType);
  potentialMaster.setPairPotential(oType, oType, &pOO);
  potentialMaster.setPairPotential(hType, hType, &pHH);
  potentialMaster.setPairPotential(hType, mType, &pHM);
  potentialMaster.setPairPotential(mType, mType, &pMM);
  potentialMaster.setCharge(hType, qH);
  potentialMaster.setCharge(mType, -2*qH);
  potentialMaster.setEwald(kCut, alpha);
  potentialMaster.setDoTruncationCorrection(false);
  potentialMaster.init();
  potentialMaster.reset();
  //PotentialMaster potentialMaster(plj, box);
  //IntegratorNVE integrator(speciesList.getAtomInfo(), potentialMaster, rand, box);
  PotentialCallbackEnergy pce;
  MeterFullCompute meterFull0(potentialMaster);
  meterFull0.setDoCompute(true);
  meterFull0.addCallback(&pce);
  for (int i=0; i<100; i++) {
    meterFull0.getData();
  }
  double* data = meterFull0.getData();
  printf("u %f %f\n", data[0], data[0]/numMolecules);
  IntegratorNHC integrator(speciesList.getAtomInfo(), potentialMaster, rand, box, 3, 0.1);
  integrator.setTimeStep(0.005);
  integrator.setTemperature(temperature);
  integrator.setNbrCheckInterval(20);
  integrator.reset();
  return 0;
  PotentialCallbackHMA pcHMA(box, temperature, 9.550752245164025e+00);
  printf("u: %f\n", integrator.getPotentialEnergy());
  if (doData) integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 1);
  MeterFullCompute meterFull(potentialMaster);
  meterFull.setDoCompute(false);
  PotentialCallbackPressure pcp(box, temperature, true);
  meterFull.addCallback(&pcp);
  integrator.addPotentialCallback(&pcp);
  if (doHMA) {
    meterFull.addCallback(&pcHMA);
    integrator.addPotentialCallback(&pcHMA);
  }
  DataPump pumpFull(meterFull, 1);
  MeterKineticEnergy meterKE;
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
  printf("time: %4.3f\n", t2-t1);
}

