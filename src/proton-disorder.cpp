/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <math.h>
#include "proton-disorder.h"
#include "box.h"
#include "action.h"
#include "random.h"
#include "alloc2d.h"

double** ProtonDisorder::go2(const char* configOname, const int num0, double* L0, const int* reps, Random& rand, const double drNbrOO, const double bondLengthOH, const double bondAngleHOH, const double offsetM) {
  SpeciesList speciesListO;
  SpeciesSimple speciesO(1, 16);
  speciesListO.add(&speciesO);
  Box boxO(speciesListO);
  boxO.setBoxSize(L0[0],L0[1],L0[2]);
  boxO.setNumMolecules(0, num0);

  ConfigurationFile configO(boxO, configOname);
  configO.go();
  const double* bs = boxO.getBoxSize();
  L0[0] = bs[0]; L0[1] = bs[1]; L0[2] = bs[2];
  Replicate::go(boxO, reps);

  return ProtonDisorder::go(boxO, rand, drNbrOO, bondLengthOH, bondAngleHOH, offsetM);
}

double** ProtonDisorder::go(Box& box, Random& rand, const double drNbrOO, const double bondLengthOH, const double bondAngleHOH, const double offsetM) {
  int N = box.getNumAtoms();
  int** pairOdCoord = (int**)malloc2D(N*2, 3, sizeof(int));
  int* nCoordsO = (int*)malloc(N*sizeof(int));
  for (int i=0; i<N; i++) nCoordsO[i] = 0;
  double dr[3];
  int nCoordTot = 0;
  double drNbr2 = drNbrOO*drNbrOO;
  // create initial assignments for protons (completely random)
  for (int i=0; i<N; i++) {
    double* ri = box.getAtomPosition(i);
    for (int j=i+1; j<N; j++) {
      double* rj = box.getAtomPosition(j);
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      box.nearestImage(dr);
      double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      if (dr2 < drNbr2) {
        if (nCoordTot==2*N) {
          fprintf(stderr, "I found too many neighboring oxygen pairs.  Please check your configuration and/or adjust the neighbor criteria.\n");
          abort();
        }
        // we found a pair of oxygens that needs a hydrogen
        int dRandCoordi = rand.nextInt(1); // does hydrogen belong to i
        nCoordsO[i] += dRandCoordi;
        nCoordsO[j] += 1-dRandCoordi;
        pairOdCoord[nCoordTot][0] = dRandCoordi ? i : j; // belongs to i|j
        pairOdCoord[nCoordTot][1] = dRandCoordi ? j : i; // this hydrogen points toward j|i
        nCoordTot++;
      }
    }
  }

  int test=0, steps=0;
  while (test != N) {
    // pick a random hydrogen
    int H = rand.nextInt(nCoordTot);
    int iO = pairOdCoord[H][0];
    int jO = pairOdCoord[H][1];
    int dCoordNew;

    // check for imbalance of hydrogens
    int dCoordOld = abs(nCoordsO[iO] - nCoordsO[jO]);
    // this H belongs to O i, compute new imbalance if we swap
    dCoordNew = abs((nCoordsO[iO] - 1) - (nCoordsO[jO] + 1));
    // swap if doing so decreases imbalance.  if imbalance unchanged, then swap 50% of the time
    if (dCoordNew < dCoordOld || (dCoordNew == dCoordOld && rand.nextInt(1) == 0)) {
      // yes, swap
      nCoordsO[iO]--;
      nCoordsO[jO]++;
      int t = pairOdCoord[H][0];
      pairOdCoord[H][0] = pairOdCoord[H][1];
      pairOdCoord[H][1] = t;
    }
    test = 0;
    steps++;
    for (int l=0; l<N; l++) {
      if (nCoordsO[l] == 2) test++;
    }
  }
  printf("Converged in %d MC steps\n", steps);
  // the two H belonging to each O
  int** OH12 = (int**)malloc2D(N, 2, sizeof(int));
  for (int i=0; i<N; i++) {
    OH12[i][0] = OH12[i][1] = -1;
  }
  for (int i=0; i<nCoordTot; i++) {
    // which O does this H belong to
    int Oi = pairOdCoord[i][0];
    int Oj = pairOdCoord[i][1];
    if (OH12[Oi][0] == -1) {
      // first H
      OH12[Oi][0] = Oj;
    }
    else {
      OH12[Oi][1] = Oj;
    }
  }

  int x = offsetM==0 ? 3 : 4;
  double** rHOM = (double**)malloc2D(N*x, 3, sizeof(double));
  double c = bondLengthOH * cos(0.5*bondAngleHOH);
  double s = bondLengthOH * sin(0.5*bondAngleHOH);
  for (int i=0; i<N; i++) {
    double* rOi = box.getAtomPosition(i);
    double bisect[3] = {0,0,0};
    for (int j=0; j<2; j++) {
      double* rOj = box.getAtomPosition(OH12[i][j]);
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = rOj[k] - rOi[k];
      box.nearestImage(dr);
      double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      double dr1 = sqrt(dr2);
      for (int k=0; k<3; k++) dr[k] /= dr1;
      // put H pointing from Oi to Oj, just need direction for now
      for (int k=0; k<3; k++) {
        rHOM[x*i+j][k] = dr[k];
        bisect[k] += dr[k];
      }
    }
    for (int k=0; k<3; k++) rHOM[x*i+2][k] = rOi[k];
    double dr2 = bisect[0]*bisect[0] + bisect[1]*bisect[1] + bisect[2]*bisect[2];
    double dr1 = sqrt(dr2);
    for (int k=0; k<3; k++) bisect[k] /= dr1;
    if (offsetM != 0) {
      // place M
      for (int k=0; k<3; k++) rHOM[x*i+3][k] = rOi[k] + bisect[k]*offsetM;
    }
    double perp[3];
    double dot = 0;
    for (int k=0; k<3; k++) dot += bisect[k]*rHOM[x*i][k];
    dr2 = 0;
    for (int k=0; k<3; k++) {
      perp[k] = rHOM[x*i][k] - bisect[k]*dot;
      dr2 += perp[k]*perp[k];
    }
    dr1 = sqrt(dr2);
    for (int k=0; k<3; k++) perp[k] /= dr1;
    for (int k=0; k<3; k++) {
      rHOM[x*i][k] = rOi[k] + bisect[k]*c + perp[k]*s;
      rHOM[x*i+1][k] = rOi[k] + bisect[k]*c - perp[k]*s;
    }
  }
  free2D((void**)pairOdCoord);
  free2D((void**)OH12);
  return rHOM;
}

void ProtonDisorder::freeRHOM(double** rHOM) {
  free2D((void**)rHOM);
}

void ProtonDisorder::snapWater(Box& box, Random& rand, const double bondLengthOH, const double bondAngleHOH, const double offsetM) {
  int N = box.getTotalNumMolecules();
  double c = bondLengthOH * cos(0.5*bondAngleHOH);
  double s = bondLengthOH * sin(0.5*bondAngleHOH);
  for (int iMolecule=0; iMolecule<N; iMolecule++) {
    int iFirstAtom = 3*iMolecule;
    double* rOi = box.getAtomPosition(iFirstAtom+2);
    double bisect[3] = {0,0,0};
    double rH0[3] = {0,0,0};
    for (int j=0; j<2; j++) {
      double* rHj = box.getAtomPosition(iFirstAtom+1+j);
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = rHj[k] - rOi[k];
      box.nearestImage(dr);
      double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      double dr1 = sqrt(dr2);
      for (int k=0; k<3; k++) dr[k] /= dr1;
      // just need direction for now
      for (int k=0; k<3; k++) {
        if (j==0) rH0[k] = dr[k];
        bisect[k] += dr[k];
      }
    }
    double dr2 = bisect[0]*bisect[0] + bisect[1]*bisect[1] + bisect[2]*bisect[2];
    double dr1 = sqrt(dr2);
    for (int k=0; k<3; k++) bisect[k] /= dr1;
    if (offsetM != 0) {
      // place M
      double* rMi = box.getAtomPosition(iFirstAtom+3);
      for (int k=0; k<3; k++) rMi[k] = rOi[k] + bisect[k]*offsetM;
    }
    double perp[3];
    double dot = 0;
    for (int k=0; k<3; k++) dot += bisect[k]*rH0[k];
    dr2 = 0;
    for (int k=0; k<3; k++) {
      perp[k] = rH0[k] - bisect[k]*dot;
      dr2 += perp[k]*perp[k];
    }
    dr1 = sqrt(dr2);
    for (int k=0; k<3; k++) perp[k] /= dr1;
    double* rH1i = box.getAtomPosition(iFirstAtom);
    double* rH2i = box.getAtomPosition(iFirstAtom+1);
    for (int k=0; k<3; k++) {
      rH1i[k] = rOi[k] + bisect[k]*c + perp[k]*s;
      rH2i[k] = rOi[k] + bisect[k]*c - perp[k]*s;
    }
  }
}
