#include "integrator.h"

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
  if (listener->callFinished) listenersStepFinished.push_back(listener);
}

void Integrator::removeListener(IntegratorListener* listener) {
  if (listener->callFinished) {
    std::vector<IntegratorListener*>::iterator it = find (listenersStepFinished.begin(), listenersStepFinished.end(), listener);
    listenersStepFinished.erase(it);
  }
}

IntegratorListener::IntegratorListener() {
  callFinished = callAccept = callReject = false;
}
