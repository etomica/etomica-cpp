#include "integrator.h"

Integrator::Integrator(PotentialMaster& p) : potentialMaster(p), temperature(1), stepCount(0) {
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
  listeners.push_back(listener);
}


IntegratorListener::IntegratorListener() {}
