/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <cstdio>

#include "potential-master.h"
#include "integrator.h"
#include "potential.h"
#include "move-virial.h"
#include "box.h"
#include "data-sink.h"
#include "random.h"
#include "util.h"
#include "cluster.h"
#include "P2PotentialGroupBuilder.h"
#include "potential-molecular.h"
#include "unit.h"
#include "virial.h"
#include "TraPPEParams.h"

int main(int argc, char** argv) {
  TraPPEParams TP;
  int nPoints = 2;
  int nDer = 0;
  double temperatureK = 1000;
  double temperature = Kelvin::toSim(1000);
  printf("temp: %f\n", temperature);
  long numSteps = 100000000;
  double sigmaHSRef = 6;
  printf("Overlap sampling for TraPPE CO2 at %.1f K for B%d and %d derivatives\n", temperatureK, nPoints, nDer);
  double vhs = 4.0/3.0*M_PI*sigmaHSRef*sigmaHSRef*sigmaHSRef;
  double HSBn = pow(vhs, nPoints-1)/2;
  for (int i=2; i<=nPoints; i++) HSBn *= i;

  Random seed;
  printf("random seed: %d\n", seed.getSeed());
  printf("Reference diagram: B%d for hard spheres with diameter %.1f Angstroms\n", nPoints, sigmaHSRef);
  printf("B%d HS: %22.15e\n", nPoints, HSBn);


  PotentialHS pHS(6);
  PotentialMolecularAtomic pHSMolecule(3, pHS);
  Box refBox(TP.speciesList);
  refBox.setBoxSize(5,5,5);
  refBox.setPeriodic(false, false, false);
  refBox.setNumMolecules(0, nPoints);
  PotentialMasterVirial potentialGroup(TP.speciesList, refBox);

  PotentialMasterVirial refPotentialMasterTraPPE = P2PotentialGroupBuilder::builder(potentialGroup, TP.speciesList, TP.numAtomTypes, TP.sigma, TP.epsilon, TP.charge, refBox);
  PotentialMasterVirialMolecular refPotentialMasterHS(TP.speciesList, refBox);
  refPotentialMasterHS.setMoleculePairPotential(0, 0, &pHSMolecule);
  IntegratorMC refIntegrator(refPotentialMasterHS, seed);
  ClusterVirial refClusterTraPPE(refPotentialMasterTraPPE, temperature, 0, false);
  ClusterChain refClusterHS(refPotentialMasterHS, temperature, 1.0, 0, false);
  MCMoveChainVirial refMove(TP.speciesList, refBox, refPotentialMasterHS, seed, sigmaHSRef);
  refIntegrator.addMove(&refMove, 1);
  MCMoveMoleculeRotateVirial refMove1(TP.speciesList, 0, refBox, refPotentialMasterHS, seed, 1.5, refClusterHS);
  // refIntegrator.addMove(&refMove1, 1);
  refIntegrator.setTemperature(temperature);

  Box targetBox(TP.speciesList);
  targetBox.setBoxSize(5,5,5);
  targetBox.setPeriodic(false, false, false);
  targetBox.setNumMolecules(0, nPoints);
  PotentialMasterVirial potentialGroupTarget(TP.speciesList, targetBox);

  PotentialMasterVirial targetPotentialMasterTraPPE = P2PotentialGroupBuilder::builder(potentialGroupTarget, TP.speciesList, TP.numAtomTypes, TP.sigma, TP.epsilon, TP.charge, targetBox);

  PotentialMasterVirialMolecular targetPotentialMasterHS(TP.speciesList, targetBox);
  targetPotentialMasterHS.setMoleculePairPotential(0, 0, &pHSMolecule);
  IntegratorMC targetIntegrator(targetPotentialMasterTraPPE, seed);
  ClusterVirial targetClusterTraPPE0(targetPotentialMasterTraPPE, temperature, 0, true);
  targetIntegrator.addListener(&targetClusterTraPPE0);
  ClusterChain targetClusterHS(targetPotentialMasterHS, temperature, 1, 0, true);
  targetIntegrator.addListener(&targetClusterHS);
  MCMoveMoleculeDisplacementVirial targetMove0(TP.speciesList, 0, targetBox, targetPotentialMasterTraPPE, seed, 1.5, targetClusterTraPPE0);
  targetIntegrator.addMove(&targetMove0, 1);
  MCMoveMoleculeRotateVirial targetMove1(TP.speciesList, 0, targetBox, targetPotentialMasterTraPPE, seed, 1.5, targetClusterTraPPE0);
  // targetIntegrator.addMove(&targetMove1, 1);
  targetIntegrator.setTemperature(temperature);

  double t1 = getTime();
  VirialAlpha *virialAlpha = new VirialAlpha(refIntegrator, targetIntegrator, refClusterHS, refClusterTraPPE, targetClusterHS, targetClusterTraPPE0);
  virialAlpha->setVerbose(true);
  virialAlpha->run();
  double t2 = getTime();

  double alpha, alphaErr, alphaCor;
  long alphaSteps = refIntegrator.getStepCount() + targetIntegrator.getStepCount();
  printf("alpha steps: %ld\n", alphaSteps);
  virialAlpha->getNewAlpha(alpha, alphaErr, alphaCor);
  printf("alpha  avg: %22.15e   err: %12.5e   cor: % 6.4f\n", alpha, alphaErr, alphaCor);
  printf("alpha time: %4.3f\n\n", t2-t1);
  // exit(0);
  long blockSize = virialAlpha->getTargetAverage().getBlockSize();
  if (blockSize > numSteps/10) {
    fprintf(stderr, "block size for uncorrelated data is large (%ld) compared to number of steps (%ld)\n", blockSize, numSteps);
  }
  delete virialAlpha;
  targetIntegrator.removeMove(&targetMove0);
  targetIntegrator.removeListener(&targetClusterTraPPE0);

  double targetStepSize = targetMove0.getStepSize();
  printf("target step size: %f\n", targetStepSize);

  ClusterVirial targetClusterTraPPE(targetPotentialMasterTraPPE, temperature, nDer, true);
  MCMoveMoleculeDisplacementVirial targetMove(TP.speciesList, 0, targetBox, targetPotentialMasterTraPPE, seed, targetStepSize, targetClusterTraPPE);
  targetIntegrator.addMove(&targetMove, 1);
  targetIntegrator.addListener(&targetClusterTraPPE);
  targetIntegrator.setTuning(false);
  VirialProduction virialProduction(refIntegrator, targetIntegrator, refClusterHS, refClusterTraPPE, targetClusterHS, targetClusterTraPPE, alpha, HSBn);
  virialProduction.runSteps(numSteps);
  double t3 = getTime();
  double acceptance = targetMove.getAcceptance();
  printf("target move acceptance: %5.3f\n", acceptance);
  virialProduction.printResults(nullptr);

  printf("time: %4.3f\n", t3-t2);
}

