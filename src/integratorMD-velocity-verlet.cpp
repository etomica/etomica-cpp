/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "integrator.h"
#include "potential-master.h"

/**
 * IntegratorMDVelocityVerlet implements velocity-verlet integration for
 * molecular dynamics.  On its own, this class samples NVE.  This class can
 * handle as NVT by adding an AndersenThermostat listener.
 */

IntegratorMDVelocityVerlet::IntegratorMDVelocityVerlet(AtomInfo& ai, PotentialMaster& p, Random& r, Box& b) : IntegratorMD(ai, p, r, b) {
}

IntegratorMDVelocityVerlet::~IntegratorMDVelocityVerlet(){}

void IntegratorMDVelocityVerlet::doStep() {
  stepCount++;
  for (vector<IntegratorListener*>::iterator it = listenersStepStarted.begin(); it!=listenersStepStarted.end(); it++) {
    (*it)->stepStarted();
  }
  int n = box.getNumAtoms();
  for (int iAtom=0; iAtom<n; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double rMass = 1.0/atomInfo.getMass(iType);
    double* ri = box.getAtomPosition(iAtom);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] += 0.5*tStep*forces[iAtom][j]*rMass;
      ri[j] += tStep*vi[j];
    }
  }

  for (vector<IntegratorListener*>::iterator it = listenersPreForce.begin(); it!=listenersPreForce.end(); it++) {
    (*it)->preForce();
  }

  computeForces();

  for (vector<IntegratorListener*>::iterator it = listenersPostForce.begin(); it!=listenersPostForce.end(); it++) {
    (*it)->postForce();
  }

  kineticEnergy = 0;
  for (int iAtom=0; iAtom<n; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double mass = atomInfo.getMass(iType);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] += 0.5*tStep*forces[iAtom][j]/mass;
      kineticEnergy += 0.5*mass*vi[j]*vi[j];
    }
  }
  for (vector<IntegratorListener*>::iterator it = listenersStepFinished.begin(); it!=listenersStepFinished.end(); it++) {
    (*it)->stepFinished();
  }
#ifdef DEBUG
  printf("%ld step %e %e %e\n", stepCount, energy/n, kineticEnergy/n, (energy + kineticEnergy)/n);
#endif
}

void IntegratorMDVelocityVerlet::reset() {
  IntegratorMD::reset();
  randomizeVelocities(false);
}
