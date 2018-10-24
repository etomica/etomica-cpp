/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "integrator.h"
#include "potential-master.h"

Integrator::Integrator(PotentialMaster& p) : potentialMaster(p), temperature(1), energy(0), stepCount(0) {
  callFinished = true;
  selfPotentialCallbackVec.push_back(this);
}

void Integrator::allComputeFinished(double uTotNew, double virialTotNew, double **f) {
  energy = uTotNew;
}

void Integrator::doSteps(int steps) {
  for (int i=0; i<steps; i++) {
    doStep();
  }
}

long Integrator::getStepCount() {
  return stepCount;
}

void Integrator::setTemperature(double T) {
  temperature = T;
}

double Integrator::getTemperature() {
  return temperature;
}

void Integrator::reset() {
  potentialMaster.computeAll(selfPotentialCallbackVec);
}

double Integrator::getPotentialEnergy() {
  return energy;
}

void Integrator::addListener(IntegratorListener* listener) {
  if (listener->callStepStarted) listenersStepStarted.push_back(listener);
  if (listener->callStepFinished) listenersStepFinished.push_back(listener);
  if (listener->callPreForce) listenersPreForce.push_back(listener);
  if (listener->callPostForce) listenersPostForce.push_back(listener);
}

void Integrator::removeListener(IntegratorListener* listener) {
  if (listener->callStepFinished) {
    std::vector<IntegratorListener*>::iterator it = find (listenersStepFinished.begin(), listenersStepFinished.end(), listener);
    listenersStepFinished.erase(it);
  }
  if (listener->callPreForce) {
    std::vector<IntegratorListener*>::iterator it = find (listenersPreForce.begin(), listenersPreForce.end(), listener);
    listenersPreForce.erase(it);
  }
  if (listener->callPostForce) {
    std::vector<IntegratorListener*>::iterator it = find (listenersPreForce.begin(), listenersPostForce.end(), listener);
    listenersPostForce.erase(it);
  }
}

IntegratorListener::IntegratorListener() {
  callStepStarted = callStepFinished = callAccept = callReject = false;
}
