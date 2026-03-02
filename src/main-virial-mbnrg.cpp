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
#include "diagrams.h"
#include "P2PotentialGroupBuilder.h"
#include "potential-angle.h"
#include "potential-molecular.h"
#include "unit.h"
#include "virial.h"
#include "TraPPEParams.h"
#include "potential/2b/energy2b.h"

int main(int argc, char** argv) {
  vector<double> xyz1 = {
    0.000000, 0.000000, 0.000000,
-0.332280, 1.099680, -0.160915,
0.332280, -1.099680, 0.160915
  };
  vector<double> xyz2 = {
    -0.062697, 1.144255, 0.511376,
-0.574729, 2.185132, 0.511376,
0.449335, 0.103379, 0.511376,
  };

  double energy = e2b::get_2b_energy("co2", "co2", 1, xyz1, xyz2);
  printf("%f \n", Energy::toSim(energy));
  TraPPEParams TP(TraPPEParams::CO2);
  int nPoints = 2;
  int nDer = 0;
  double temperatureK = 1000;
  double temperature = Kelvin::toSim(1000);
  printf("temp: %f\n", temperature);
  long numSteps = 100000000;
  double sigmaHSRef = 6;
  printf("Overlap sampling for TraPPE %d at %.1f K for B%d and %d derivatives\n", TP.chemForm,temperatureK, nPoints, nDer);
  double vhs = 4.0/3.0*M_PI*sigmaHSRef*sigmaHSRef*sigmaHSRef;
  double HSBn = pow(vhs, nPoints-1)/2;
  for (int i=2; i<=nPoints; i++) HSBn *= i;
  bool anyPolar = false, anyFlex = false;
  anyPolar = anyPolar || TP.isPolar;
  anyFlex = anyFlex || TP.isFlex;

  Random seed;
  printf("random seed: %d\n", seed.getSeed());
  printf("Reference diagram: B%d for hard spheres with diameter %.1f Angstroms\n", nPoints, sigmaHSRef);
  printf("B%d HS: %22.15e\n", nPoints, HSBn);


  PotentialHS pHS(6);
  PotentialMolecularAtomic pHSMolecule(TP.speciesList, pHS);
  Box refBox(TP.speciesList);
  refBox.setBoxSize(5,5,5);
  refBox.setPeriodic(false, false, false);
  refBox.setNumMolecules(0, nPoints);
  // PotentialMasterVirial potentialGroup(TP.speciesList, refBox);
  //flex moves

  PotentialMasterIntraMBnrg potentialMasterIntraRef(refBox);

  int nSpheres = TP.species.getNumAtoms();
  printf("nSpheres = %d\n", nSpheres);
  bool isFlex = TP.isFlex && (TP.diagram.empty() || TP.diagram != "BC");
  printf("isFlex = %s\n", isFlex ? "true" : "false");
  printf("Diagram %s\n" , TP.diagram.c_str());

  // PotentialMasterVirial refPotentialMasterTraPPE = P2PotentialGroupBuilder::builder(potentialGroup, TP.speciesList, TP.species.getNumAtomTypes(), TP.sigma, TP.epsilon, TP.charge, refBox);
  PotentialMolecularMBnrg potentialGroup;
  std::vector<double> xyz1ut = {6.6630444410e-01,  -3.8357176030e-01, 1.1519802350e-01,  1.7838183644e+00,
                                    -1.9222069500e-01, -1.7587628680e-01, -4.4811475090e-01, -5.7997649630e-01,
                                    4.1069507510e-01};
  std::vector<double> xyz2ut = {2.4803292099e+00,  7.5103875900e-01,  -2.9043390394e+00, 2.2674715176e+00,
                                    1.8651909097e+00,  -2.6082983571e+00, 2.7020706245e+00,  -3.5351972400e-01,
                                    -3.2106052081e+00};
  for (int i = 0; i < 3; i++)
  {
    double* r1 = refBox.getAtomPosition(i);
    double* r2 = refBox.getAtomPosition(3+i);

    for (int j = 0; j < 3; j++)
    {
      r1[j] = xyz1ut[3*i+j];
      r2[j] = xyz2ut[3*i+j];
    }

  }
  double e = e2b::get_2b_energy("co2_archive", "co2_archive", 1, xyz1ut, xyz2ut);
  printf("%f \n", e);
  potentialGroup.u(refBox, 0, 1);
  exit(0);
  PotentialMasterVirialMolecular refPotentialMasterTraPPE(TP.speciesList, refBox);
  refPotentialMasterTraPPE.setMoleculePairPotential(0, 0, &potentialGroup);


  PotentialMasterVirialMolecular refPotentialMasterHS(TP.speciesList, refBox);
  refPotentialMasterHS.setMoleculePairPotential(0, 0, &pHSMolecule);
  IntegratorMC refIntegrator(refPotentialMasterHS, seed);
  ClusterVirial refClusterTraPPE(refPotentialMasterTraPPE, temperature, 0, false);
  ClusterChain refClusterHS(refPotentialMasterHS, temperature, 1.0, 0, false);
  MCMoveChainVirial refMove(TP.speciesList, refBox, seed, sigmaHSRef);
  refIntegrator.addMove(&refMove, 1);
  MCMoveMoleculeRotateVirial refMove1(TP.speciesList, 0, refBox, seed, 1.5, refClusterHS);
  refIntegrator.addMove(&refMove1, 1);
  refIntegrator.setTemperature(temperature);
  MCMoveClusterAngleGeneral refMoveAngle(refBox, potentialMasterIntraRef, false, seed, TP.bonding, TP.triplets, TP.speciesList, 0.1, refClusterHS);
  refIntegrator.addMove(&refMoveAngle, 1);
  refMoveAngle.idx = 0;
  MCMoveClusterStretch refMoveStretch(refBox, potentialMasterIntraRef, seed, TP.bonding, TP.speciesList, 0.1, refClusterHS);
  refIntegrator.addMove(&refMoveStretch, 1);

  Box targetBox(TP.speciesList);
  targetBox.setBoxSize(5,5,5);
  targetBox.setPeriodic(false, false, false);
  targetBox.setNumMolecules(0, nPoints);
  PotentialMolecularMBnrg potentialGroupTarget;
  PotentialMasterVirialMolecular targetPotentialMasterTraPPE(TP.speciesList, targetBox);
  targetPotentialMasterTraPPE.setMoleculePairPotential(0, 0, &potentialGroupTarget);

  //flex moves

  PotentialMasterIntraMBnrg potentialMasterIntraTarget(targetBox);

  //P2 Stretch

  potentialMasterIntraRef.init();
  potentialMasterIntraTarget.init();
  PotentialMasterVirialMolecular targetPotentialMasterHS(TP.speciesList, targetBox);
  targetPotentialMasterHS.setMoleculePairPotential(0, 0, &pHSMolecule);
  IntegratorMC targetIntegrator(targetPotentialMasterTraPPE, seed);
  Cluster* targetClusterTraPPE0 = nullptr;
  if (anyFlex && nPoints==2 && anyPolar)
  {
    targetClusterTraPPE0 = new ClusterVirial(targetPotentialMasterTraPPE, temperature, 0, false);

    //flipping for flexible polar B2
    vector<vector <int>> flipPoints = diagrams::getFlipPointsforDiagram("1");
    targetClusterTraPPE0 = new ClusterFlippedPoints(*targetClusterTraPPE0, TP.speciesList, targetBox, false, flipPoints, 3);
  }
  else
  {
    targetClusterTraPPE0 = new ClusterVirial(targetPotentialMasterTraPPE, temperature, 0, true);
  }

  targetIntegrator.addListener(targetClusterTraPPE0);
  ClusterChain targetClusterHS(targetPotentialMasterHS, temperature, 1, 0, true);
  targetIntegrator.addListener(&targetClusterHS);
  MCMoveMoleculeDisplacementVirial targetMove0(TP.speciesList, 0, targetBox, seed, 1.5, *targetClusterTraPPE0);
  targetIntegrator.addMove(&targetMove0, 1);
  MCMoveMoleculeRotateVirial targetMove1(TP.speciesList, 0, targetBox, seed, 1.5, *targetClusterTraPPE0);
  targetIntegrator.addMove(&targetMove1, 1);
  MCMoveClusterAngleGeneral targetMoveAngle(targetBox, potentialMasterIntraTarget, false, seed, TP.bonding, TP.triplets, TP.speciesList, 0.1, *targetClusterTraPPE0);
  targetIntegrator.addMove(&targetMoveAngle, 1);
  targetMoveAngle.idx = 1;
  MCMoveClusterStretch targetMoveStretch(targetBox, potentialMasterIntraTarget, seed, TP.bonding, TP.speciesList, 0.1, *targetClusterTraPPE0);
  targetIntegrator.addMove(&targetMoveStretch, 1);
  targetMoveStretch.idx = 1;

  // for (int i=3; i<=5 ; i++)
  // {
  //   targetBox.getAtomPosition(i)[1]+=5.5;
  // }
  targetIntegrator.setTemperature(temperature);
  // for (int i = 0 ; i <= 10 ; i ++)
  // {
  //   targetIntegrator.doStep();
  // }
  // exit(0);
  double t1 = getTime();
  VirialAlpha *virialAlpha = new VirialAlpha(refIntegrator, targetIntegrator, refClusterHS, refClusterTraPPE, targetClusterHS, *targetClusterTraPPE0);
  virialAlpha->setVerbose(true);
  virialAlpha->run();
  double t2 = getTime();

  double alpha, alphaErr, alphaCor;
  long alphaSteps = refIntegrator.getStepCount() + targetIntegrator.getStepCount();
  printf("alpha steps: %ld\n", alphaSteps);
  virialAlpha->getNewAlpha(alpha, alphaErr, alphaCor);
  printf("alpha  avg: %22.15e   err: %12.5e   cor: % 6.4f\n", alpha, alphaErr, alphaCor);
  printf("alpha time: %4.3f\n\n", t2-t1);
  // alpha = 1;
  // exit(0);
  long blockSize = virialAlpha->getTargetAverage().getBlockSize();
  if (blockSize > numSteps/10) {
    fprintf(stderr, "block size for uncorrelated data is large (%ld) compared to number of steps (%ld)\n", blockSize, numSteps);
  }
  delete virialAlpha;
  targetIntegrator.removeMove(&targetMove0);
  targetIntegrator.removeMove(&targetMove1);
  targetIntegrator.removeMove(&targetMoveAngle);
  targetIntegrator.removeMove(&targetMoveStretch);

  targetIntegrator.removeListener(targetClusterTraPPE0);

  double targetStepSize = targetMove0.getStepSize();
  printf("target step size: %f\n", targetStepSize);
  Cluster* targetClusterTraPPE = nullptr;
  if (anyFlex && nPoints==2 && anyPolar)
  {
    targetClusterTraPPE = new ClusterVirial(targetPotentialMasterTraPPE, temperature, nDer, false);

    //flipping for flexible polar B2
    vector<vector <int>> flipPoints = diagrams::getFlipPointsforDiagram("1");
    targetClusterTraPPE = new ClusterFlippedPoints(*targetClusterTraPPE, TP.speciesList, targetBox, false, flipPoints, 3);
  }
  else
  {
    targetClusterTraPPE = new ClusterVirial(targetPotentialMasterTraPPE, temperature, nDer, true);

  }

  MCMoveMoleculeDisplacementVirial targetMove(TP.speciesList, 0, targetBox, seed, targetStepSize, *targetClusterTraPPE);
  MCMoveMoleculeRotateVirial targetRotateMove(TP.speciesList, 0, targetBox, seed, targetStepSize, *targetClusterTraPPE);
  MCMoveClusterAngleGeneral targetAngleMove(targetBox, potentialMasterIntraTarget, false, seed, TP.bonding, TP.triplets, TP.speciesList, 0.1, *targetClusterTraPPE);
  MCMoveClusterStretch targetStretchMove(targetBox, potentialMasterIntraTarget, seed, TP.bonding, TP.speciesList, 0.1, *targetClusterTraPPE);
  targetStretchMove.idx = 1;
  targetIntegrator.addMove(&targetRotateMove, 1);
  // targetAngleMove.idx = 2;
  targetIntegrator.addMove(&targetMove, 1);
  targetIntegrator.addMove(&targetAngleMove, 1);
  targetIntegrator.addMove(&targetStretchMove, 1);

  targetIntegrator.addListener(targetClusterTraPPE);
  targetIntegrator.setTuning(false);
  VirialProduction virialProduction(refIntegrator, targetIntegrator, refClusterHS, refClusterTraPPE, targetClusterHS, *targetClusterTraPPE, alpha, HSBn);
  virialProduction.runSteps(numSteps);
  double t3 = getTime();
  double acceptance = targetMove.getAcceptance();
  printf("target move acceptance: %5.3f\n", acceptance);
  virialProduction.printResults(nullptr);
  // for (int i=0; i<=5 ; i++)
  // {
  //   double *r = targetBox.getAtomPosition(i);
  //   printf("%f %f %f \n", r[0], r[1], r[2]);
  // }

  printf("time: %4.3f\n", t3-t2);
}

