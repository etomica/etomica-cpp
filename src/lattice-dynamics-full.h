/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifdef LATTICE_DYNAMICS
#pragma once

#include <vector>
#include "potential-callback.h"

using namespace std;

class PotentialMaster;
class PotentialCallback;

class LatticeDynamicsFull : public PotentialCallbackMoleculePhi {
  protected:
    vector<PotentialCallback*> selfPotentialCallbackVec;
    int numCells[3];
    int wCount;
    double** waveVectors;
    std::complex<double> ***matrix;
    double** evals;
    int* wvCount;
    double logSum;
    bool unstable;
  public:
    LatticeDynamicsFull(PotentialMaster& potentialMaster);
    virtual ~LatticeDynamicsFull();
    virtual void allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor);
    void setNumCells(int x, int y, int z);
    void setupForWV(int n, double** wv);
    void compute();
    int getNumWavevectors() {return wCount;}
    double** getWavevectors() {return waveVectors;}
    bool getUnstable() {return unstable;}
    double getLogSum() {return logSum;}
    double** getEigenvalues() {return evals;}
};
#endif
