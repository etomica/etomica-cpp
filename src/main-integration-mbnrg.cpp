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
  double temperatureK = 400;
  long step = 10;
  double minDist = 0;
  double maxDist = 7;
  double stepSize = 6;
  long intraSteps = 10;
  bool rigid = false;

  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if (arg == "--temp" && i + 1 < argc) {
      temperatureK = std::stod(argv[++i]);
    } else if (arg == "--step" && i + 1 < argc) {
      step = std::stol(argv[++i]);
    } else if (arg == "--minDist" && i + 1 < argc) {
      minDist = std::stod(argv[++i]);
    } else if (arg == "--maxDist" && i + 1 < argc) {
      maxDist = std::stod(argv[++i]);
    } else if (arg == "--intraSteps" && i + 1 < argc) {
      intraSteps = std::stol(argv[++i]);
    } else if (arg == "--stepSize" && i + 1 < argc) {
      stepSize = std::stod(argv[++i]);
    }
    else if (arg == "--rigid" && i + 1 < argc) {
      rigid = std::stoi(argv[++i]);
    }




  }
  double temperature = Kelvin::toSim(temperatureK);
  std::string compound_name = TraPPEParams::chemFormToString(TP.chemForm);
  printf("Integration for MBnrg %s at %.1f K for B2 with dist: %f-%f and step size: %f\n", compound_name.c_str(),temperatureK, minDist, maxDist, stepSize);
  std::cout << "Steps: " << step << "\n";
  std::cout << "Intra Steps: " << intraSteps << "\n";
  std::cout << "rigid: " << rigid << "\n";

  bool anyFlex = TP.isFlex && !rigid;

  Random seed = 1639844264;
  printf("random seed: %d\n", seed.getSeed());

  bool isFlex = TP.isFlex && (TP.diagram.empty() || TP.diagram != "BC");
  printf("isFlex = %s\n", isFlex ? "true" : "false");

  PotentialMolecularMBnrg potentialGroup;
  PotentialHS pHS(1000);
  PotentialMolecularAtomic pHSMolecule(TP.speciesList, pHS);

  // exit(0);
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

  potentialMasterIntraTarget.init();
  PotentialMasterVirialMolecular targetPotentialMasterHS(TP.speciesList, targetBox);
  targetPotentialMasterHS.setMoleculePairPotential(0, 0, &pHSMolecule);
  IntegratorMC targetIntegrator(targetPotentialMasterTraPPE, seed);
  Cluster* targetClusterTraPPE0 = nullptr;
  PotentialMolecularMBnrg3body p3BodyTarget;
  ClusterChain dummyCluster(targetPotentialMasterHS, temperature, 1, 0, false);
  if (anyFlex && nPoints==2)
  {
    targetClusterTraPPE0 = new ClusterMBX(targetPotentialMasterTraPPE, temperature, 0, false, p3BodyTarget);

    //flipping for flexible polar B2
    vector<vector <int>> flipPoints = diagrams::getFlipPointsforDiagram("1");
    targetClusterTraPPE0 = new ClusterFlippedPoints(*targetClusterTraPPE0, TP.speciesList, targetBox, false, flipPoints, 3);
  }
  else
  {
    //B3 biconnected
    targetClusterTraPPE0 = new ClusterMBX(targetPotentialMasterTraPPE, temperature, 0, false, p3BodyTarget);
    targetClusterTraPPE0 = new ClusterFlipped(*targetClusterTraPPE0, TP.speciesList, targetBox, false);
  }

  targetIntegrator.addListener(targetClusterTraPPE0);
  MCMoveMoleculeRotateVirial targetMove1(TP.speciesList, 0, targetBox, seed, 0.1, dummyCluster);
  targetIntegrator.addMove(&targetMove1, 0.1);
  MCMoveClusterAngleGeneral targetMoveAngle(targetBox, potentialMasterIntraTarget, false, seed, TP.bonding, TP.triplets, TP.speciesList, 0.1, dummyCluster);
  if (anyFlex) targetIntegrator.addMove(&targetMoveAngle, 1);
  targetMoveAngle.idx = 1;
  MCMoveClusterStretch targetMoveStretch(targetBox, potentialMasterIntraTarget, seed, TP.bonding, TP.speciesList, 0.1, dummyCluster);
  if (anyFlex) targetIntegrator.addMove(&targetMoveStretch, 1);
  targetMoveStretch.idx = 1;

  for (int i=3; i<=5 ; i++)
  {
    targetBox.getAtomPosition(i)[0];
  }
  targetIntegrator.setTemperature(temperature);
  // for (int i = 0 ; i <= 10 ; i ++)
  // {
  //   targetIntegrator.doStep();
  // }
  // exit(0);
  int steps = round((maxDist-minDist)/stepSize)+1;
  double* sums = new double[steps];
  double* SS = new double[steps];
  double totalIntegral = 0;
  double SStotalIntegral = 0;
  for (int l=0; l<steps; l++) sums[l]=SS[l]=0;
  for (int i=0; i < step; i++) {
    for (int j=0; j < intraSteps; j++) {
      targetIntegrator.doStep();
    }
    double thisIntegral = 0;
    for (int k=0;k<steps;k++) {
      double x2 = minDist + k*stepSize;
      double disp2[3] = {x2 - targetBox.getAtomPosition(3)[0],-targetBox.getAtomPosition(3)[1] , -targetBox.getAtomPosition(3)[2]};
      for (int a=3; a<6; a++ ) for (int m=0; m<3; m++)targetBox.getAtomPosition(a)[m]+=disp2[m];

      double disp1[3] = {-targetBox.getAtomPosition(0)[0],-targetBox.getAtomPosition(0)[1] , -targetBox.getAtomPosition(0)[2]};
      for (int a=0; a<3; a++ )for (int m=0; m<3; m++)targetBox.getAtomPosition(a)[m]+=disp1[m];
      double value = targetClusterTraPPE0->getValues()[0]; //integrand
      sums[k]+=value;
      double fac = (k==0||k==steps-1)?0.5:1;
      thisIntegral +=fac*value*4*M_PI*x2*x2*stepSize;
      SS[k]+=value*value;


    }
    totalIntegral +=thisIntegral;
    SStotalIntegral += thisIntegral*thisIntegral;
    printf("%d steps have finished\n", i+1);
  }

  for (int k=0;k<steps;k++) printf("%f %f %f\n", minDist + k*stepSize, sums[k]/step, 0);
  double avg = totalIntegral/step;
  double err = 0;
  printf("Integral = %f %f\n", avg, err);
}