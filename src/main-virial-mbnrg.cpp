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
int main(int argc, char** argv) {
  TraPPEParams TP(TraPPEParams::CO2);
  int nPoints = 2;
  int nDer = 0;
  double temperatureK = 400;
  long numSteps = 10000;
  double discreteStepSize = 6;
  double discreteCutOff = 7;
  bool fixedCOMflag = false;
  bool rigid = true;
  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if (arg == "--temp" && i + 1 < argc) {
      temperatureK = std::stod(argv[++i]);
    } else if (arg == "--steps" && i + 1 < argc) {
      numSteps = std::stol(argv[++i]);
    } else if (arg == "--nPoints" && i + 1 < argc) {
      nPoints = std::stoi(argv[++i]);
    }
    else if (arg == "--discreteStepSize" && i + 1 < argc) {
    discreteStepSize = std::stod(argv[++i]);
    }
    else if (arg == "--discreteCutOff" && i + 1 < argc) {
      discreteCutOff = std::stod(argv[++i]);
    }
    else if (arg == "--fixedCOMflag" && i + 1 < argc) {
      fixedCOMflag = std::stoi(argv[++i]);
    }
    else if (arg == "--rigid" && i + 1 < argc) {
      rigid = std::stoi(argv[++i]);
    }



  }
  double temperature = Kelvin::toSim(temperatureK);
  double sigmaHSRef = 6.025;
  std::string compound_name = TraPPEParams::chemFormToString(TP.chemForm);
  printf("Overlap sampling for MBnrg %s at %.1f K for B%d and %d derivatives\n", compound_name.c_str(),temperatureK, nPoints, nDer);
  std::cout << "Steps: " << numSteps << "\n";
  std::cout << "Discrete CutOff: " << discreteCutOff << "\n";
  std::cout << "Discrete StepSize: " << discreteStepSize << "\n";
  std::cout << "FixedCOM flag: " << fixedCOMflag << "\n";
  std::cout << "rigid: " << rigid << "\n";


  double vhs = 4.0/3.0*M_PI*sigmaHSRef*sigmaHSRef*sigmaHSRef;
  printf("vhs: %f\n", vhs);
  if (discreteStepSize > 0) {
    vhs = 0;
    int n = sigmaHSRef/discreteStepSize;
    int pSum = n * (n+1) * (1+2*n) / 6;

    vhs = 4.0 * M_PI * discreteStepSize * discreteStepSize * discreteStepSize * pSum;
    printf("discrete vhs: %f\n", vhs);

  }
  // exit(0);
  double HSBn = pow(vhs, nPoints-1)/2;


  for (int i=2; i<=nPoints; i++) HSBn *= i;
  bool anyPolar = TP.isPolar;
  bool anyFlex = TP.isFlex && !rigid;

  Random seed = 1639844264;
  printf("random seed: %d\n", seed.getSeed());
  printf("Reference diagram: B%d for hard spheres with diameter %.3f Angstroms\n", nPoints, sigmaHSRef);
  printf("B%d HS: %22.15e\n", nPoints, HSBn);


  PotentialHS pHS(sigmaHSRef);
  PotentialMolecular* pHSMolecule;
  if (discreteStepSize > 0) pHSMolecule = new PotentialMolecularAtomicC(TP.speciesList, pHS);
  else
  {
    pHSMolecule = new PotentialMolecularAtomic(TP.speciesList, pHS);
  }
  Box refBox(TP.speciesList);
  refBox.idx = 0;
  refBox.setBoxSize(5,5,5);
  refBox.setPeriodic(false, false, false);
  refBox.setNumMolecules(0, nPoints);
  // PotentialMasterVirial potentialGroup(TP.speciesList, refBox);
  //flex moves

  PotentialMasterIntraMBnrg potentialMasterIntraRef(refBox);
  int nSpheres = TP.species.getNumAtoms();
  printf("nSpheres = %d\n", nSpheres);
  //flexible corrections
  bool isFlex = TP.isFlex && (TP.diagram.empty() || TP.diagram != "BC");
  printf("isFlex = %s\n", isFlex ? "true" : "false");
  printf("Diagram %s\n" , TP.diagram.c_str());

  // PotentialMasterVirial refPotentialMasterTraPPE = P2PotentialGroupBuilder::builder(potentialGroup, TP.speciesList, TP.species.getNumAtomTypes(), TP.sigma, TP.epsilon, TP.charge, refBox);
  PotentialMolecularMBnrg potentialGroup;
  // std::vector<double> xyz1ut = {6.6630444410e-01,  -3.8357176030e-01, 1.1519802350e-01,  1.7838183644e+00,
  //                                   -1.9222069500e-01, -1.7587628680e-01, -4.4811475090e-01, -5.7997649630e-01,
  //                                   4.1069507510e-01};
  // std::vector<double> xyz2ut = {2.4803292099e+00,  7.5103875900e-01,  -2.9043390394e+00, 2.2674715176e+00,
  //                                   1.8651909097e+00,  -2.6082983571e+00, 2.7020706245e+00,  -3.5351972400e-01,
  //                                   -3.2106052081e+00};
  // for (int i = 0; i < 3; i++)
  // {
  //   double* r1 = refBox.getAtomPosition(i);
  //   double* r2 = refBox.getAtomPosition(3+i);
  //
  //   for (int j = 0; j < 3; j++)
  //   {
  //     r1[j] = xyz1ut[3*i+j];
  //     r2[j] = xyz2ut[3*i+j];
  //   }
  //
  // }
  // exit(0);
  PotentialMasterVirialMolecular refPotentialMasterTraPPE(TP.speciesList, refBox);
  refPotentialMasterTraPPE.setMoleculePairPotential(0, 0, &potentialGroup);


  PotentialMasterVirialMolecular refPotentialMasterHS(TP.speciesList, refBox);
  refPotentialMasterHS.setMoleculePairPotential(0, 0, pHSMolecule);
  IntegratorMC refIntegrator(refPotentialMasterHS, seed);
  PotentialMolecularMBnrg3body p3BodyRef;
  ClusterMBX refClusterTraPPE(refPotentialMasterTraPPE, temperature, 0, false, p3BodyRef);
  ClusterChain refClusterHS(refPotentialMasterHS, temperature, 1.0, 0, false);
  MCMoveChainVirial refMove(TP.speciesList, refBox, seed, sigmaHSRef, discreteStepSize);
  refIntegrator.addMove(&refMove, 1);
  refMove.fixedCOM = fixedCOMflag;
  MCMoveMoleculeRotateVirial refMove1(TP.speciesList, 0, refBox, seed, 1.5, refClusterHS);
  refIntegrator.addMove(&refMove1, 1);
  refIntegrator.setTemperature(temperature);
  refMove1.fixedCOM = fixedCOMflag;
  MCMoveClusterAngleGeneral refMoveAngle(refBox, potentialMasterIntraRef, false, seed, TP.bonding, TP.triplets, TP.speciesList, 0.1, refClusterHS);
  if (anyFlex) refIntegrator.addMove(&refMoveAngle, 1);
  refMoveAngle.fixedCOM = fixedCOMflag;
  refMoveAngle.idx = 0;
  MCMoveClusterStretch refMoveStretch(refBox, potentialMasterIntraRef, seed, TP.bonding, TP.speciesList, 0.1, refClusterHS);
  if (anyFlex) refIntegrator.addMove(&refMoveStretch, 1);
  refMoveStretch.fixedCOM = fixedCOMflag;

  Box targetBox(TP.speciesList);
  targetBox.idx = 1;
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
  targetPotentialMasterHS.setMoleculePairPotential(0, 0, pHSMolecule);
  IntegratorMC targetIntegrator(targetPotentialMasterTraPPE, seed);
  Cluster* targetClusterTraPPE0 = nullptr;
  PotentialMolecularMBnrg3body p3BodyTarget;
  printf("anyPolar %d \n", anyPolar);
  if (nPoints==2 && anyPolar)
  {
    targetClusterTraPPE0 = new ClusterMBX(targetPotentialMasterTraPPE, temperature, 0, false, p3BodyTarget);

    //flipping for flexible polar B2
    vector<vector <int>> flipPoints = diagrams::getFlipPointsforDiagram("1");
    if (!rigid) targetClusterTraPPE0 = new ClusterFlippedPoints(*targetClusterTraPPE0, TP.speciesList, targetBox, false, flipPoints, 3);
  }
  else
  {
    //B3 biconnected
    targetClusterTraPPE0 = new ClusterMBX(targetPotentialMasterTraPPE, temperature, 0, false, p3BodyTarget);
    targetClusterTraPPE0 = new ClusterFlipped(*targetClusterTraPPE0, TP.speciesList, targetBox, true);
  }

  targetIntegrator.addListener(targetClusterTraPPE0);
  ClusterChain targetClusterHS(targetPotentialMasterHS, temperature, 1, 0, true);
  targetIntegrator.addListener(&targetClusterHS);
  MCMoveMoleculeDisplacementVirial targetMove0(TP.speciesList, 0, targetBox, seed, 6.0, *targetClusterTraPPE0, discreteCutOff, discreteStepSize);

  targetMove0.tunable = false;
  targetIntegrator.addMove(&targetMove0, 1);
  MCMoveMoleculeRotateVirial targetMove1(TP.speciesList, 0, targetBox, seed, 1.5, *targetClusterTraPPE0);
  targetIntegrator.addMove(&targetMove1, 1);
  targetMove1.fixedCOM = fixedCOMflag;
  MCMoveClusterAngleGeneral targetMoveAngle(targetBox, potentialMasterIntraTarget, false, seed, TP.bonding, TP.triplets, TP.speciesList, 0.1, *targetClusterTraPPE0);
  if (anyFlex) targetIntegrator.addMove(&targetMoveAngle, 1);
  targetMoveAngle.idx = 1;
  targetMoveAngle.fixedCOM = fixedCOMflag;
  MCMoveClusterStretch targetMoveStretch(targetBox, potentialMasterIntraTarget, seed, TP.bonding, TP.speciesList, 0.1, *targetClusterTraPPE0);
  if (anyFlex) targetIntegrator.addMove(&targetMoveStretch, 1);
  targetMoveStretch.idx = 1;
  targetMoveStretch.fixedCOM = fixedCOMflag;

  // for (int i=3; i<=5 ; i++)
  // {
  //   targetBox.getAtomPosition(i)[0];
  // }

  targetIntegrator.setTemperature(temperature);
  // for (int i=3; i<=5 ; i++)
  // {
  //   targetBox.getAtomPosition(i)[1]+=5.5;
  // }

  // for (int i = 0 ; i <= 10 ; i ++)
  // {
  //   targetIntegrator.doStep();
  // }
  // exit(0);

  double t1 = getTime();

  VirialAlpha *virialAlpha = new VirialAlpha(refIntegrator, targetIntegrator, refClusterHS, refClusterTraPPE, targetClusterHS, *targetClusterTraPPE0);
  virialAlpha->targetMeter.box = &targetBox;
  virialAlpha->refMeter.box = &refBox;

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
  if (anyFlex) targetIntegrator.removeMove(&targetMoveAngle);
  if (anyFlex) targetIntegrator.removeMove(&targetMoveStretch);

  targetIntegrator.removeListener(targetClusterTraPPE0);

  double targetStepSize = targetMove0.getStepSize();
  printf("target step size: %f\n", targetStepSize);
  Cluster* targetClusterTraPPE = nullptr;
  printf("anyPOlar=%d \n", anyPolar);
  if (nPoints==2 && anyPolar)
  {
    targetClusterTraPPE = new ClusterMBX(targetPotentialMasterTraPPE, temperature, nDer, false, p3BodyTarget);

    //flipping for flexible polar B2
    vector<vector <int>> flipPoints = diagrams::getFlipPointsforDiagram("1");
    if (!rigid) targetClusterTraPPE = new ClusterFlippedPoints(*targetClusterTraPPE, TP.speciesList, targetBox, false, flipPoints, 3);
  }
  else
  {
    targetClusterTraPPE = new ClusterMBX(targetPotentialMasterTraPPE, temperature, nDer, true, p3BodyTarget);

  }

  MCMoveMoleculeDisplacementVirial targetMove(TP.speciesList, 0, targetBox, seed, targetStepSize, *targetClusterTraPPE, discreteCutOff, discreteStepSize);
  MCMoveMoleculeRotateVirial targetRotateMove(TP.speciesList, 0, targetBox, seed, targetStepSize, *targetClusterTraPPE);
  targetRotateMove.fixedCOM = fixedCOMflag;
  MCMoveClusterAngleGeneral targetAngleMove(targetBox, potentialMasterIntraTarget, false, seed, TP.bonding, TP.triplets, TP.speciesList, 0.1, *targetClusterTraPPE);
  targetAngleMove.fixedCOM = fixedCOMflag;
  MCMoveClusterStretch targetStretchMove(targetBox, potentialMasterIntraTarget, seed, TP.bonding, TP.speciesList, 0.1, *targetClusterTraPPE);
  targetStretchMove.idx = 1;
  targetStretchMove.fixedCOM = fixedCOMflag;
  targetIntegrator.addMove(&targetRotateMove, 1);
  // targetAngleMove.idx = 2;
  targetIntegrator.addMove(&targetMove, 1);
  if (anyFlex) targetIntegrator.addMove(&targetAngleMove, 1);
  if (anyFlex) targetIntegrator.addMove(&targetStretchMove, 1);

  targetIntegrator.addListener(targetClusterTraPPE);
  targetIntegrator.setTuning(false);
  MeterVirialDirect directMeter(refClusterHS, refClusterTraPPE);
  Average directAverage(1, 10, 1000, false);
  DataPump directPump(directMeter, 1, &directAverage);
  refIntegrator.addListener(&directPump);

  VirialProduction virialProduction(refIntegrator, targetIntegrator, refClusterHS, refClusterTraPPE, targetClusterHS, *targetClusterTraPPE, alpha, HSBn);
  virialProduction.runSteps(numSteps);
  printf("%f %f\n", directAverage.getStatistics()[0][AVG_AVG], directAverage.getStatistics()[0][AVG_ERR]);
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

