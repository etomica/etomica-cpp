/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <math.h>
#include "proton-disorder.h"
#include "box.h"
#include "random.h"
#include "alloc2d.h"

double** ProtonDisorder::go(Box& box, Random& rand, const double drNbrOO, const double bondLengthOH, const double bondAngleHOH, const double offsetM) {
  int N = box.getNumAtoms();
  int** pairOdCoord = (int**)malloc2D(N*2, 3, sizeof(int));
  int* nCoordsO = (int*)malloc(N*sizeof(int));
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

  int x = offsetM==0 ? 2 : 3;
  double** rH = (double**)malloc2D(N*x, 3, sizeof(double));
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
        rH[x*i+j][k] = dr[k];
        bisect[k] += dr[k];
      }
    }
    double dr2 = bisect[0]*bisect[0] + bisect[1]*bisect[1] + bisect[2]*bisect[2];
    double dr1 = sqrt(dr2);
    for (int k=0; k<3; k++) bisect[k] /= dr1;
    if (offsetM != 0) {
      // place M
      for (int k=0; k<3; k++) rH[x*i+2][k] = rOi[k] + bisect[k]*offsetM;
    }
    double perp[3];
    double dot = 0;
    for (int k=0; k<3; k++) dot += bisect[k]*rH[x*i][k];
    dr2 = 0;
    for (int k=0; k<3; k++) {
      perp[k] = rH[x*i][k] - bisect[k]*dot;
      dr2 += perp[k]*perp[k];
    }
    dr1 = sqrt(dr2);
    for (int k=0; k<3; k++) perp[k] /= dr1;
    for (int k=0; k<3; k++) {
      rH[x*i][k] = rOi[k] + bisect[k]*c + perp[k]*s;
      rH[x*i+1][k] = rOi[k] + bisect[k]*c - perp[k]*s;
    }
  }
  return rH;
}
