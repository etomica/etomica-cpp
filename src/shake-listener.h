/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "integrator.h"
class RigidConstraint;

class ShakeListener : public IntegratorListener {
  protected:
    SpeciesList& speciesList;
    Box& box;
    IntegratorMD& integrator;
    vector<RigidConstraint*> *constraints;
    double tol;
    int maxIterations;
    int maxNumAtoms, maxNumBonds;
    double **rOld, **rNew, **rijOld, **rijNew;
    double *dl, *dd, *du, *dl0, *dd0, *du0, *rhs, *rhsOld;
    double *lambda;
    void tridiagSolve(int n, double* dl, double* dd, double* du, double* rhs);

  public:
    ShakeListener(SpeciesList& speciesList, Box& box, IntegratorMD& integrator);
    virtual ~ShakeListener();
    virtual void init();
    virtual void stepStarted();
    virtual void preForce();
    virtual void stepFinished();
};
