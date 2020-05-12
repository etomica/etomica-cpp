/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifdef LATTICE_DYNAMICS
#pragma once

class Potential;

class LatticeDynamics {
  private:
    int numCells[3];
    double density;
    double rCut;
    double lrcFac;
    Potential* potential;
    bool doLRC;
    int status;
    int nBasis;
    double **basis;
    double cellSize[3];
    double strain[3];
    double ***selfSum;
    int wCount;
    double** waveVectors;
    std::complex<double> ***matrix;
    int* wvCount;
    int latticeShells[3];
    long long doneXYZ;

    double rNext;
    int xcellNext, ycellNext, zcellNext;
    int wvNext;

    double logSum;
    bool unstable;
    double uLat;
  public:
    LatticeDynamics(double d, Potential *p, bool doLRC, int nBasis);
    ~LatticeDynamics();
    void setBasis(int i, double x, double y, double z);
    void setUnitCell(double x, double y, double z);
    void setStrain(double x, double y, double z);
    void setNumCells(int x, int y, int z);
    int getStatus() {return status;}
    void setup();
    void setupForWV(int n, double** wv);
    long long countLS();
    long long goLS(int nMax);
    int doSelfSum();
    int goEVD(int nwv);
    // returned eigenvalues should be freed by caller
    double** getEigenvalues(int nwv, int& rwv);
    int getWaveVectorCount() {return wCount;}
    bool getUnstable() {return unstable;}
    double getLogSum() {return logSum;}
    double getU();
#endif
};
