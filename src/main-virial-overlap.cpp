#include <stdio.h>

#include "potential-master.h"
#include "integrator.h"
#include "potential.h"
#include "move-virial.h"
#include "box.h"
#include "meter-virial.h"
#include "data-sink.h"
#include "random.h"
#include "util.h"
#include "cluster.h"
#include "virial.h"

int main(int argc, char** argv) {
  int order = 4;
  double temperature = 1.0;
  long steps = 10000000;
  double sigmaRef = 1.5;

  double vSphere = 4.0/3.0*M_PI*sigmaRef*sigmaRef*sigmaRef;
  double refIntegral = pow(vSphere, order-1)/2;
  for (int i=2; i<=order; i++) refIntegral *= i;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());
  printf("Reference Integral: %22.15e\n", refIntegral);

  PotentialLJ plj(1,1,TRUNC_NONE, 1.0/0.0);
  PotentialHS pHS(1.5);
  SpeciesList speciesList;
  speciesList.add(new SpeciesSimple(1,1));

  Box refBox(speciesList);
  refBox.setBoxSize(5,5,5);
  refBox.setNumMolecules(0, order);
  PotentialMasterVirial refPotentialMasterLJ(speciesList, refBox);
  refPotentialMasterLJ.setPairPotential(0, 0, &plj);
  PotentialMasterVirial refPotentialMasterHS(speciesList, refBox);
  refPotentialMasterHS.setPairPotential(0, 0, &pHS);
  IntegratorMC refIntegrator(refPotentialMasterHS, rand);
  ClusterVirial refClusterLJ(refPotentialMasterLJ, temperature, 0, false);
  ClusterChain refClusterHS(refPotentialMasterHS, temperature, 1.0, 0);
  MCMoveChainVirial refMove(refBox, refPotentialMasterHS, rand, 1.5);
  refIntegrator.addMove(&refMove, 1);
  refIntegrator.setTemperature(temperature);

  Box targetBox(speciesList);
  targetBox.setBoxSize(5,5,5);
  targetBox.setNumMolecules(0, order);
  PotentialMasterVirial targetPotentialMasterLJ(speciesList, targetBox);
  targetPotentialMasterLJ.setPairPotential(0, 0, &plj);
  PotentialMasterVirial targetPotentialMasterHS(speciesList, targetBox);
  targetPotentialMasterHS.setPairPotential(0, 0, &pHS);
  IntegratorMC targetIntegrator(targetPotentialMasterLJ, rand);
  ClusterVirial targetClusterLJ(targetPotentialMasterLJ, temperature, 0, true);
  ClusterChain targetClusterHS(targetPotentialMasterHS, temperature, 1, 0);
  MCMoveDisplacementVirial targetMove(targetBox, targetPotentialMasterLJ, rand, 0.2, targetClusterLJ);
  targetIntegrator.addMove(&targetMove, 1);
  targetIntegrator.setTemperature(temperature);

  double t1 = getTime();
  VirialAlpha *virialAlpha = new VirialAlpha(refIntegrator, targetIntegrator, refClusterHS, refClusterLJ, targetClusterHS, targetClusterLJ);
  virialAlpha->run();
  double t2 = getTime();

  double alpha, alphaErr, alphaCor;
  long alphaSteps = refIntegrator.getStepCount() + targetIntegrator.getStepCount();
  printf("alpha steps: %ld\n", alphaSteps);
  virialAlpha->getNewAlpha(alpha, alphaErr, alphaCor);
  printf("alpha  avg: %22.15e   err: %12.5e   cor: % 6.4f\n", alpha, alphaErr, alphaCor);
  printf("alpha time: %4.3f\n\n", t2-t1);
  long blockSize = virialAlpha->getTargetAverage().getBlockSize();
  if (blockSize > steps/10) {
    fprintf(stderr, "block size for uncorrelated data is large (%ld) compared to number of steps (%ld)\n", blockSize, steps);
  }
  delete virialAlpha;

  targetIntegrator.setTuning(false);
  refIntegrator.setTuning(false);
  double targetStepSize = targetMove.getStepSize();
  printf("target step size: %f\n", targetStepSize);

  VirialProduction virialProduction(refIntegrator, targetIntegrator, refClusterHS, refClusterLJ, targetClusterHS, targetClusterLJ, alpha, refIntegral);
  virialProduction.runSteps(steps, steps/1000);
  double t3 = getTime();
  double acceptance = targetMove.getAcceptance();
  printf("target move acceptance: %5.3f\n", acceptance);
  const char *targetNames[] = {"B"};
  virialProduction.printResults(targetNames);

  printf("time: %4.3f\n", t3-t2);
}

