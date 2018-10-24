/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "integrator.h"

IntegratorMulti::IntegratorMulti(PotentialMaster& pm) : Integrator(pm) {
  callFinished = false;
}

IntegratorMulti::~IntegratorMulti() {}

void IntegratorMulti::addIntegrator(Integrator* i, long iSteps) {
  integrators.push_back(i);
  steps.push_back(iSteps);
}

void IntegratorMulti::setSteps(Integrator* i, long iSteps) {
  for (int j=0; j<(int)integrators.size(); j++) {
    if (integrators[j] == i) {
      steps[j] = iSteps;
      return;
    }
  }
  fprintf(stderr, "Can't find the integrator you specified\n");
  abort();
}

void IntegratorMulti::doStep() {
  stepCount++;
  for (vector<IntegratorListener*>::iterator it = listenersStepStarted.begin(); it!=listenersStepStarted.end(); it++) {
    (*it)->stepStarted();
  }
  for (int j=0; j<(int)integrators.size(); j++) {
    integrators[j]->reset();
    integrators[j]->doSteps(steps[j]);
  }
  for (vector<IntegratorListener*>::iterator it = listenersStepFinished.begin(); it!=listenersStepFinished.end(); it++) {
    (*it)->stepFinished();
  }
}

void IntegratorMulti::reset() {
  for (int j=0; j<(int)integrators.size(); j++) {
    integrators[j]->reset();
  }
}

Integrator* IntegratorMulti::getIntegrator(int j) {
  if (j<0 || j>(int)integrators.size()) {
    fprintf(stderr, "Invalid integrator %d\n", j);
    abort();
  }
  return integrators[j];
}
