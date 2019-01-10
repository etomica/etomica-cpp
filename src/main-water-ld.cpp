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
#include "util.h"
#include "minimize.h"
#include "action.h"
#include "ewald.h"
#include "proton-disorder.h"
#include "lattice-dynamics-full.h"
#include "parameter-map.h"
#include "alloc2d.h"

void vecE(Box& box, int iAtom, double* v) {
  double* r = box.getAtomPosition(iAtom);
  r[0] = v[0];
  r[1] = v[1];
  r[2] = v[2];
  box.nearestImage(r);
}

int main(int argc, char** argv) {
  ParameterMap pm;
  pm.addParameter("N0", "46");
  pm.addParameter("L0", "12.03");
  pm.addParameter("disorderRep", "1");
  pm.addParameter("ldRep", "1");
  pm.addParameter("uTol", "5e-8");
  pm.addParameter("s", "4");
  pm.addParameter("alpha", "0.3");
  pm.addParameter("steps", "20");
  if (argc>1) pm.parseArgs(argc-1, argv+1);
  int N0 = pm.getInt("N0");
  double L0 = pm.getDouble("L0");
  double density = 0;
  if (pm.hasParameter("density")) {
    density = pm.getDouble("density");
    L0 = pow(N0/density, 1.0/3.0);
  }
  int disorderRep = pm.getInt("disorderRep");
  int ldRep = pm.getInt("ldRep");
  bool cheap = pm.hasParameter("cheap") && pm.getBool("cheap");
  double s = pm.getDouble("s");

  Random rand;
  if (pm.hasParameter("seed")) {
    rand.setSeed(pm.getInt("seed"));
  }
  printf("random seed: %d\n", rand.getSeed());

  int numMolecules = N0*disorderRep*disorderRep*disorderRep;
  double uTol = pm.getDouble("uTol");
  double L0k[3] = {L0, L0, L0};
  int rep[3] = {disorderRep, disorderRep, disorderRep};
  double** rHOM = ProtonDisorder::go2("O.pos", N0, L0k, rep, rand, 3.5, 0.9572, 104.52/180*M_PI, 0.15);

  SpeciesList speciesList;
  SpeciesFile species("water.species");
  speciesList.add(&species);

  double L[3] = {L0k[0]*disorderRep, L0k[1]*disorderRep, L0k[2]*disorderRep};
  printf("box size: %f %f %f\n", L[0], L[1], L[2]);
  Box box(speciesList);
  box.setNumMolecules(0, numMolecules);
  box.setBoxSize(L[0],L[1],L[2]);
  for (int i=0; i<numMolecules; i++) {
    vecE(box, 4*i, rHOM[4*i+0]);
    vecE(box, 4*i+1, rHOM[4*i+1]);
    vecE(box, 4*i+2, rHOM[4*i+2]);
    vecE(box, 4*i+3, rHOM[4*i+3]);
  }

  if (density > 0) {
    double scaleL = pow(numMolecules/density / (L[0]*L[1]*L[2]), 1.0/3.0);
    L[0] *= scaleL; L[1] *= scaleL; L[2] *= scaleL;
    box.scaleBoxTo(L[0],L[1],L[2]);
  }

  double A = 600E3, C = 610;
  double s6 = A/C;
  double sigma = pow(s6, 1.0/6.0);
  double epsilon = (C/s6)*1000/4 * 1000.*1e20*1e-24*4.184;
  EwaldFourier ewald(speciesList, box);
  double alpha, rc, kCut;
  if (cheap) {
    if (!pm.hasParameter("rCut")) {
      fprintf(stderr, "rCut must be specified if using cheap potential\n");
      abort();
    }
    rc = pm.getDouble("rCut");
    alpha = pm.getDouble("alpha");
    kCut = 0;
  }
  else {
    ewald.getOptimalAlpha(s, alpha, rc, kCut);
    double f = pow(numMolecules*4, 1.0/6.0);
    alpha /= f;
    rc *= f;
    kCut /= f;
  }
  printf("alpha: %f\n", alpha);
  printf("rc: %f\n", rc);
  printf("kc: %f\n", kCut);
  PotentialSS pOO12(4*epsilon*pow(sigma,12),12,TRUNC_SIMPLE,rc);
  double qH = 193.82504408037946;
  PotentialEwaldBare pHH(alpha, qH*qH, rc, TRUNC_FORCE_SHIFT);
  PotentialEwaldBare pHM(alpha, -2*qH*qH, rc, TRUNC_FORCE_SHIFT);
  PotentialEwaldBare pMM(alpha, 4*qH*qH, rc, TRUNC_FORCE_SHIFT);
  double uq = 4*qH*qH*erfc(rc*alpha)/rc;
  double eta = PotentialEwald6::getEta(rc, sigma, epsilon, uq);
  printf("eta: %f\n", eta);
  PotentialEwald6 pOO(pOO12, sigma, epsilon, sigma, epsilon, eta, rc, TRUNC_FORCE_SHIFT);
  int oType = species.getTypeForSymbol("O");
  int mType = species.getTypeForSymbol("M");
  int hType = species.getTypeForSymbol("H");

  int nc = 3;
  for (int k=0; k<3; k++) {
    if (rc/L[k] > nc) nc = ((int)rc/L[k]) + 1;
  }
  PotentialMasterCell potentialMaster(speciesList, box, false, nc);
  potentialMaster.setPairPotential(oType, oType, &pOO);
  potentialMaster.setPairPotential(hType, hType, &pHH);
  potentialMaster.setPairPotential(hType, mType, &pHM);
  potentialMaster.setPairPotential(mType, mType, &pMM);
  ewald.setCharge(hType, qH);
  ewald.setCharge(mType, -2*qH);
  ewald.setCutoff(kCut);
  ewald.setChargeAlpha(alpha);
  ewald.setR6Coeff(oType, sigma, epsilon);
  ewald.setR6eta(eta);
  potentialMaster.setEwald(&ewald);
  potentialMaster.setDoTruncationCorrection(false);
  potentialMaster.init();
  int* numCells = potentialMaster.getNumCells();
  printf("cells(%d): %d %d %d\n", nc, numCells[0], numCells[1], numCells[2]);
  PotentialCallbackEnergy pce;
  MeterFullCompute meterFull(potentialMaster);
  meterFull.addCallback(&pce);
  double* data = meterFull.getData();
  double lastU = data[0]/numMolecules;
  printf("u0: %f\n", lastU);

  double t1 = getTime();
  Minimize min(potentialMaster, false);
  int steps = pm.getInt("steps");
  char** symbols = (char**)malloc2D(3,2,sizeof(char));
  symbols[0][0] = 'H';
  symbols[1][0] = 'O';
  symbols[2][0] = 'M';
  symbols[0][1] = symbols[1][1] = symbols[2][1] = '\0';
  if (pm.hasParameter("xyz") && pm.getBool("xyz")) {
    WriteXYZ::go("water.xyz", box, symbols, false);
  }
  for (int i=0; i<steps; i++) {
    min.doStep();
    if (pm.hasParameter("xyz") && pm.getBool("xyz")) {
      WriteXYZ::go("water.xyz", box, symbols, true);
    }
    double lastDR = min.getLastDR();
    data = meterFull.getData();
    double unow = data[0]/numMolecules;
    double du = unow-lastU;
    printf("%d % 10.4e % 10.4e\n", i+1, lastDR, unow-lastU);
    lastU = unow;
    if (fabs(du) < uTol) break;
  }
  double t2 = getTime();
  printf("u: %20.15e\n", lastU);
  PotentialCallbackPressure pcp(box, 0, true);
  meterFull.addCallback(&pcp);
  PotentialCallbackVirialTensor pcvt(box);
  meterFull.addCallback(&pcvt);
  data = meterFull.getData();
  printf("p: %20.15e\n", data[1]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", data[2], data[3], data[4]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", data[3], data[5], data[6]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", data[4], data[6], data[7]);
  printf("min time: %4.3f\n", t2-t1);

  printf("=== Lattice dynamics ===\n");
  if (ldRep>1) {
    int ldReplicates[3] = {ldRep,ldRep,ldRep};
    Replicate::go(box, ldReplicates);
    printf("box size: %f\n", box.getBoxSize()[0]);
    numMolecules *= ldRep*ldRep*ldRep;

    if (!cheap) {
      alpha /= ldRep;
      rc *= ldRep;
      kCut /= ldRep*ldRep;
      pOO.setCutoff(rc);
      pHH.setCutoff(rc);
      pHM.setCutoff(rc);
      pMM.setCutoff(rc);
      ewald.setCutoff(kCut);
      ewald.setChargeAlpha(alpha);
      uq = 4*qH*qH*erfc(rc*alpha)/rc;
      eta = PotentialEwald6::getEta(rc, sigma, epsilon, uq);
      printf("alpha: %f\n", alpha);
      printf("rc: %f\n", rc);
      printf("kc: %f\n", kCut);
      printf("eta: %f\n", eta);
    }
  }
  potentialMaster.init();
  LatticeDynamicsFull ld(potentialMaster);
  ld.setNumCells(ldRep,ldRep,ldRep);
  t1 = getTime();
  ld.compute();
  t2 = getTime();

  if  (ld.getUnstable()) {
    printf("System is unstable\n");
    exit(1);
  }
  double logSum = ld.getLogSum();
  printf("lnsum: %f\n", logSum/numMolecules);
  printf("ld time: %f\n", t2-t1);
}

