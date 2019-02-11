/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <vector>
#include "potential-callback.h"

using namespace std;

class PotentialMaster;
class PotentialCallback;

class Minimize : public PotentialCallbackMoleculePhi {
  protected:
    long stepCount;
    vector<PotentialCallback*> selfPotentialCallbackVec;
    double* fMolecule;
    double lastDR;
    double lastU, lastDU;
    double lastStep;
    double maxDR, maxDtheta, maxDL;
    bool sdStep;
  public:
    Minimize(PotentialMaster& potentialMaster, bool flexBox);
    virtual ~Minimize();
    void doStep();
    void doSteps(int steps, double uTol, double rTol);
    long getStepCount();
    virtual void allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor);
    double getLastDR();
    double getLastU();
    double getLastDU();
};

class OptimizeBoxShape : public PotentialCallback {
  protected:
    PotentialMaster& potentialMaster;
    long stepCount;
    vector<PotentialCallback*> selfPotentialCallbackVec;
    double lastU, lastDU;
    double stepSize;
    double virialTensor[6];
  public:
    OptimizeBoxShape(PotentialMaster& potentialMaster);
    virtual ~OptimizeBoxShape();
    void doStep();
    void doSteps(int steps, double uTol);
    long getStepCount();
    virtual void allComputeFinished(double uTot, double virialTot, double** f, double* virialTensor);
    double getLastU();
    double getLastDU();
};
