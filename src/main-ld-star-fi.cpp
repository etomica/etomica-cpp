#include <stdio.h>
#include <Eigen/Dense>
#include "ld.h"
#include "potential.h"
#include "alloc2d.h"

void usage() {
  fprintf(stderr, "usage: ld-star-fi-star [--ss|--ss6] [--lrc|--shift|--fshift] [--atomicExp] ncells density rc nwv\n");
}

int main(int argc, char** argv) {
  if (argc < 5 || argc > 8) {
    usage();
    exit(1);
  }
  int arg0 = 1;
  bool ss = false, ss6 = false;
  if (!strcmp(argv[1], "--ss")) {
    arg0++;
    ss = true;
  }
  else if (!strcmp(argv[1], "--ss6")) {
    arg0++;
    ss6 = true;
  }
  int trunc = TRUNC_SIMPLE;
  bool lrc = false;
  if (!strcmp(argv[arg0], "--lrc")) {
    lrc = true;
    arg0++;
  }
  else if (!strcmp(argv[arg0], "--shift")) {
    trunc = TRUNC_SHIFT;
    printf(" shifted\n");
    arg0++;
  }
  else if (!strcmp(argv[arg0], "--fshift")) {
    trunc = TRUNC_FORCE_SHIFT;
    printf(" force-shifted\n");
    arg0++;
  }

  bool atomicExp = false;
  if (!strcmp(argv[arg0], "--atomicExp")) {
    arg0++;
    atomicExp = true;
    printf(" doing atomic exponential: exp(-i k.r_ij)\n");
  }

  int nCells = atoi(argv[arg0]);
  char* c2 = argv[arg0+1];
  char* c1 = c2;
  double density = strtod(c1, &c2);
  if (c1==c2) {
    usage();
    exit(1);
  }
  c2 = argv[arg0+2];
  c1 = c2;
  double rc = strtod(c1, &c2);
  if (c1==c2) {
    usage();
    exit(1);
  }
  int nwv = atoi(argv[arg0+3]);
  Potential* p = nullptr;
  if (ss) p = new PotentialSS(4.0, 12, trunc, rc);
  else if (ss6) p = new PotentialSS(4.0, 6, trunc, rc);
  else p = new PotentialLJ(1.0, 1.0, trunc, rc);
  int nBasis = 4;
  LatticeDynamics* ld = new LatticeDynamics(density, p, lrc, nBasis, atomicExp);
  ld->setNumCells(nCells,nCells,nCells);
  ld->setBasis(0, 0, 0, 0);
  ld->setBasis(1, 0.5, 0.5, 0);
  ld->setBasis(2, 0.5, 0, 0.5);
  ld->setBasis(3, 0, 0.5, 0.5);

  double cellSize = pow(nBasis/density, 1.0/3.0);
  // wvMax = kMax * wvBasis
  // kMax = nCells / 2
  // wvBasis = 2 pi / (nCells * cellSize)
  // wvMax = pi / cellSize

  double wvMax = 2.0*M_PI / cellSize;
  long nwv_tot = pow(nwv, 3);
  double** waveVectors = (double**)malloc2D(nwv_tot, 3, sizeof(double));
  int l = 0;
  for (int i=0; i<nwv; i++) {
    for (int j=0; j<nwv; j++) {
      for (int k=0; k<nwv; k++) {
        waveVectors[l][0] = i*wvMax/nwv;
        waveVectors[l][1] = j*wvMax/nwv;
        waveVectors[l][2] = k*wvMax/nwv;
        l++;
      }
    }
  }

  ld->setupForWV(nwv_tot, waveVectors);
  int n = 10000;
  while (1) {
    long long next = ld->goLS(n);
    if (next<0) break;
  }
  ld->doSelfSum();
  n=1000;
  while (1) {
    int next = ld->goEVD(n);
    if (next<0) break;
  }
  if (ld->getUnstable()) {
    printf(" unstable!\n");
  }
  else {
    long N_eff = nwv_tot*nBasis;
    double A = (ld->getLogSum())/2.0;
    double u0 = ld->getU();
    printf(" N_eff  %ld \n", N_eff);
    printf(" nc  %d \n", nCells);
    printf(" a  %15.10f \n", cellSize);
    printf(" u0 %22.15f\n", u0);
    printf(" A/N_eff  %22.15f\n", A/N_eff);
    printf(" A/N_eff-ln(N_eff)/N_eff  %22.15f\n", A/N_eff - log(N_eff)/N_eff);
  }
  delete ld;
  delete p;
}
