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
    bool sdStep;
  public:
    Minimize(PotentialMaster& potentialMaster);
    virtual ~Minimize();
    void doStep();
    void doSteps(int steps, double uTol, double rTol);
    long getStepCount();
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    double getLastDR();
    double getLastU();
    double getLastDU();
};
