/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
#include "ewald.h"

int main(int argc, char** argv) {
  int numMolecules = 46;
  double K = 0.8314459861448581;
  double temperatureK = 120;
  if (argc>1) temperatureK = atof(argv[1]);
  double temperature = temperatureK*K;
  long steps = 5000000;
  bool doData = true;
  bool doHMA = true;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  double A = 600E3, C = 610;
  double s6 = A/C;
  double sigma = pow(s6, 1.0/6.0);
  double epsilon = (C/s6)*1000/4 * 1000.*1e20*1e-24*4.184;
  printf("temperature(K): %f\n", temperatureK);
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
  EwaldFourier ewald(speciesList, box);
  ewald.setCharge(hType, qH);
  ewald.setCharge(mType, -2*qH);
  ewald.setCutoff(kCut);
  ewald.setChargeAlpha(alpha);
  potentialMaster.setEwald(&ewald);
  potentialMaster.setDoTruncationCorrection(false);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);
  PotentialCallbackEnergy pce;
  PotentialCallbackPressure pcp(box, temperature, true);
  IntegratorMC integrator(potentialMaster, rand);
  integrator.setTemperature(temperature);
  integrator.reset();
  MCMoveMoleculeDisplacement move(box, potentialMaster, rand, 0.2);
  integrator.addMove(&move, 1);
  MCMoveMoleculeRotate moveRotate(speciesList, box, potentialMaster, rand);
  integrator.addMove(&moveRotate, 1);
  double Pharm = 18.3 * (temperature/(120*K));
  PotentialCallbackMoleculeHMA pcHMA(box, speciesList, &potentialMaster, temperature, Pharm);
  pcHMA.setReturnAnharmonic(true);
  pcHMA.findShiftV();
  if (false) {
  double** dRdV = pcHMA.getDRDV();
  printf("dRdV\n");
  for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
    printf("dRdV %2d  % 12.8f % 12.8f % 12.8f  % 12.8f % 12.8f % 12.8f\n", iMolecule, dRdV[iMolecule][0], dRdV[iMolecule][1], dRdV[iMolecule][2],
                                                         dRdV[iMolecule][3], dRdV[iMolecule][4], dRdV[iMolecule][5]);
  }
  /*double* center = species.getMoleculeCOM(box, 0, 3);
  double theta = 0.0001, dx = 0.0001;
  for (int iAtom=0; iAtom<4; iAtom++) {
    double* ri = box.getAtomPosition(iAtom);
    if (true) {
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = ri[k] - center[k];
      box.nearestImage(dr);
      double dr1 = cos(theta)*dr[1] - sin(theta)*dr[2];
      double dr2 = cos(theta)*dr[2] + sin(theta)*dr[1];
      ri[1] = center[1] + dr1;
      ri[2] = center[2] + dr2;
    }
    else {
      ri[0] += dx;
    }
    box.nearestImage(ri);
  }
  potentialMaster.init();
  
  PotentialCallbackMoleculeHMA pcHMA3(box, speciesList, &potentialMaster, temperature, Pharm);
  pcHMA3.setReturnAnharmonic(true);
  pcHMA3.findShiftV();
  exit(1);*/
  int na = box.getNumAtoms();
  double L0 = L;
  L += 0.0001;
  /*for (int i=0; i<na; i++) {
    double* ri = box.getAtomPosition(i);
    for (int k=0; k<3; k++) ri[k] = ri[k]*L/L0;
  }*/
  for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
    double* center = species.getMoleculeCOM(box, iMolecule*4, iMolecule*4+3);
    for (int iAtom=iMolecule*4; iAtom<(iMolecule+1)*4; iAtom++) {
      double* ri = box.getAtomPosition(iAtom);
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = ri[k] - center[k];
      box.nearestImage(dr);
      for (int k=0; k<3; k++) ri[k] = center[k]*L/L0 + dr[k];
      //for (int k=0; k<3; k++) ri[k] = ri[k]*L/L0;
      //for (int k=0; k<3; k++) ri[k] -= dr[k]*(L-L0)/L0;
    }
  }
  box.setBoxSize(L,L,L);
  for (int i=0; i<na; i++) {
    box.nearestImage(box.getAtomPosition(i));
  }
  printf("volume scaled\n");
  printf("=================\n");
  potentialMaster.init();
  PotentialCallbackMoleculeHMA pcHMA2(box, speciesList, &potentialMaster, temperature, Pharm);
  pcHMA2.setReturnAnharmonic(true);
  pcHMA2.findShiftV();
  double dV = L*L*L-L0*L0*L0;
  printf("dV %e\n", dV);
  for (int iMolecule=0; iMolecule<numMolecules; iMolecule++) {
    double* idRdV = dRdV[iMolecule];
    double a2 = 0;
    for (int k=0; k<3; k++) a2 += idRdV[3+k]*idRdV[3+k];
    double a1 = sqrt(a2);
    double theta = asin(a1)*dV;
    double axis[3];
    for (int k=0; k<3; k++) axis[k] = a1==0 ? (k==0?1:0) : idRdV[3+k]/a1;
    double* center = species.getMoleculeCOM(box, iMolecule*4, iMolecule*4+3);
    RotationMatrix rot;
    printf("%d  % f % f % f  % 12.10f\n", iMolecule, axis[0], axis[1], axis[2], theta);
    //if (iMolecule==0) theta = 0.00000089;
    //axis[0] = 0.02; axis[1] = -0.095; axis[2] = -sqrt(1-axis[1]*axis[1]-axis[0]*axis[0]);
    rot.setAxisAngle(axis, theta);
    for (int iAtom=iMolecule*4; iAtom<(iMolecule+1)*4; iAtom++) {
      double* ri = box.getAtomPosition(iAtom);
      //printf("%d  %f %f %f  %f %f %f ", iAtom, ri[0], ri[1], ri[2], center[0], center[1], center[2]);
      rot.transformAbout(ri, center, box);
      //printf(" %f %f %f ", ri[0], ri[1], ri[2]);
      for (int k=0; k<3; k++) ri[k] += idRdV[k]*dV;
      box.nearestImage(ri);
      //printf(" %f %f %f\n", ri[0], ri[1], ri[2]);
    }
  }
  printf("lattice sites shifted\n");
  potentialMaster.init();
  PotentialCallbackMoleculeHMA pcHMA3(box, speciesList, &potentialMaster, temperature, Pharm);
  pcHMA3.setReturnAnharmonic(true);
  pcHMA3.findShiftV();
    exit(1);
  }
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
  else {
    meterFull.addCallback(&pcp);
  }
  DataPump pumpFull(meterFull, 4*numMolecules);
  if (doData) {
    integrator.addListener(&pumpPE);
    integrator.addListener(&pumpFull);
  }

  double t1 = getTime();
  integrator.doSteps(steps);
  double t2 = getTime();
  if (doData) {
    double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
    statsPE[AVG_AVG] /= numMolecules;
    statsPE[AVG_ERR] /= numMolecules;
    printf("u avg: %f  err: %f  cor: %f\n", statsPE[AVG_AVG], statsPE[AVG_ERR], statsPE[AVG_ACOR]);
    double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();
    if (!doHMA) {
      double* statsP = statsFull[0];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);
    }
    if (doHMA) {
      printf("HMA\n");
      double* statsU = statsFull[0];
      printf("u avg: %f  err: %f  cor: %f\n", statsU[AVG_AVG]/numMolecules, statsU[AVG_ERR]/numMolecules, statsU[AVG_ACOR]);
      double* statsP = statsFull[1];
      printf("p avg: %f  err: %f  cor: %f\n", statsP[AVG_AVG], statsP[AVG_ERR], statsP[AVG_ACOR]);

      double* statsUHMA = statsFull[2];
      statsUHMA[AVG_AVG] /= numMolecules;
      statsUHMA[AVG_ERR] /= numMolecules;
      printf("uHMA avg: %f  err: %f  cor: %f\n", statsUHMA[AVG_AVG], statsUHMA[AVG_ERR], statsUHMA[AVG_ACOR]);
      double* statsPHMA = statsFull[3];
      printf("pHMA avg: %f  err: %f  cor: %f\n", statsPHMA[AVG_AVG], statsPHMA[AVG_ERR], statsPHMA[AVG_ACOR]);

      double* statsFoo = statsFull[4];
      printf("pHMA1 avg: %f  err: %f  cor: %f\n", statsFoo[AVG_AVG], statsFoo[AVG_ERR], statsFoo[AVG_ACOR]);

    }
  }
  printf("time: %4.3f\n", t2-t1);
}

