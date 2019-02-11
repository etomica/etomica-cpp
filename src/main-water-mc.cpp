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
#include "data-pump.h"

void vecE(Box& box, int iAtom, double* v) {
  double* r = box.getAtomPosition(iAtom);
  r[0] = v[0];
  r[1] = v[1];
  r[2] = v[2];
  box.centralImage(r);
}

int main(int argc, char** argv) {
  ParameterMap pm;
  pm.addParameter("N0", "46");
  pm.addParameter("L0", "12.03");
  pm.addParameter("disorderRep", "1");
  pm.addParameter("mcRep", "1");
  pm.addParameter("uTol", "1e-7");
  pm.addParameter("temperatureK", "120");
  pm.addParameter("steps", "10000");
  pm.addParameter("s", "4");
  pm.addParameter("alpha", "0.3");
  if (argc>1) pm.parseArgs(argc-1, argv+1);
  int N0 = pm.getInt("N0");
  double L0 = pm.getDouble("L0");
  double density = 0;
  if (pm.hasParameter("density")) {
    density = pm.getDouble("density");
    L0 = pow(N0/density, 1.0/3.0);
  }
  int disorderRep = pm.getInt("disorderRep");
  int mcRep = pm.getInt("mcRep");
  bool cheap = pm.hasParameter("cheap") && pm.getBool("cheap");
  double s = pm.getDouble("s");
  double K = 0.8314459861448581;
  double temperatureK = pm.getDouble("temperatureK");
  double temperature = temperatureK * K;
  long steps = pm.getLong("steps");

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
  printf("cells (%d): %d %d %d\n", nc, numCells[0], numCells[1], numCells[2]);
  PotentialCallbackEnergy pce;
  MeterFullCompute meterFullE(potentialMaster);
  meterFullE.addCallback(&pce);
  double* data = meterFullE.getData();
  double lastU = data[0]/numMolecules;
  printf("u0: %f\n", lastU);
  MeterFullCompute meterFoo(potentialMaster);
  PotentialCallbackPressure pcp(box, temperature, true);
  meterFoo.addCallback(&pcp);
  PotentialCallbackVirialTensor pcvt(box);
  meterFoo.addCallback(&pcvt);
  double *p = meterFoo.getData();
  printf("p0: %f\n", p[0]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", p[1], p[2], p[3]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", p[2], p[4], p[5]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", p[3], p[5], p[6]);

  double t1 = getTime();
  Minimize min(potentialMaster, true);
  const double* bs = box.getBoxSize();
  for (int i=0; i<500; i++) {
    min.doStep();
    double lastDR = min.getLastDR();
    data = meterFullE.getData();
    double unow = data[0]/numMolecules;
    double du = unow-lastU;
    printf("%d % 10.4e % 10.4e  %8.5f %8.5f %8.5f\n", i, lastDR, unow-lastU, bs[0], bs[1], bs[2]);
    lastU = unow;
    //if (fabs(du) < uTol) break;
  }
  double t2 = getTime();
  printf("u: %20.15e\n", lastU);
  printf("min time: %4.3f\n", t2-t1);

  printf("=== Monte Carlo ===\n");
  printf("temperature %f K\n", temperatureK);
  if (mcRep>1) {
    int mcReplicates[3] = {mcRep,mcRep,mcRep};
    Replicate::go(box, mcReplicates);
  }
  if (!cheap) {
    // we optimized for 2nd derivative.  now optimize for energy
    ewald.getOptimalAlpha(s, alpha, rc, kCut);
    uq = 4*qH*qH*erfc(rc*alpha)/rc;
    eta = PotentialEwald6::getEta(rc, sigma, epsilon, uq);
    ewald.setCutoff(kCut);
    ewald.setChargeAlpha(alpha);
    ewald.setR6eta(eta);
    printf("box size: %f\n", box.getBoxSize()[0]);
    printf("alpha: %f\n", alpha);
    printf("rc: %f\n", rc);
    printf("kc: %f\n", kCut);
    printf("eta: %f\n", eta);
  }

  PotentialEwald6 pOOL(pOO12, sigma, epsilon, sigma, epsilon, eta, rc, TRUNC_FORCE_SHIFT);
  PotentialEwaldBare pHHL(alpha, qH*qH, rc, TRUNC_FORCE_SHIFT);
  PotentialEwaldBare pHML(alpha, -2*qH*qH, rc, TRUNC_FORCE_SHIFT);
  PotentialEwaldBare pMML(alpha, 4*qH*qH, rc, TRUNC_FORCE_SHIFT);
  potentialMaster.setPairPotential(oType, oType, &pOOL);
  potentialMaster.setPairPotential(hType, hType, &pHHL);
  potentialMaster.setPairPotential(hType, mType, &pHML);
  potentialMaster.setPairPotential(mType, mType, &pMML);
  potentialMaster.init();

  t1 = getTime();

  numCells = potentialMaster.getNumCells();
  printf("cells: %d %d %d\n", numCells[0], numCells[1], numCells[2]);
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

  p = meterFoo.getData();
  printf("p0: %f\n", p[0]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", p[1], p[2], p[3]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", p[2], p[4], p[5]);
  printf("  % 15.9e  % 15.9e  % 15.9e\n", p[3], p[5], p[6]);
  double u0 = integrator.getPotentialEnergy();
  printf("u: %f\n", u0/numMolecules);
  integrator.doSteps(steps/10);
  MeterPotentialEnergy meterPE(integrator);
  DataPump pumpPE(meterPE, 1);
  MeterFullCompute meterFull(potentialMaster);
  meterFull.setDoCompute(true);
  meterFull.addCallback(&pcHMA);
  DataPump pumpFull(meterFull, 4*numMolecules);

  integrator.addListener(&pumpPE);
  integrator.addListener(&pumpFull);

  printf("%ld steps\n", steps);
  t1 = getTime();
  integrator.doSteps(steps);
  t2 = getTime();

  double* statsPE = ((Average*)pumpPE.getDataSink(0))->getStatistics()[0];
  statsPE[AVG_AVG] /= numMolecules;
  statsPE[AVG_ERR] /= numMolecules;
  printf("u avg: %f  err: %f  cor: %f\n", statsPE[AVG_AVG], statsPE[AVG_ERR], statsPE[AVG_ACOR]);
  double** statsFull = ((Average*)pumpFull.getDataSink(0))->getStatistics();

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

  t2 = getTime();

  printf("mc time: %f\n", t2-t1);
}

