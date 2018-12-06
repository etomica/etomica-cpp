/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "integrator.h"
#include "potential-master.h"
#include "neighbor-update-listener.h"

IntegratorMD::IntegratorMD(AtomInfo& ai, PotentialMaster& p, Random& r, Box& b) : Integrator(p), atomInfo(ai), random(r), box(b), forces(nullptr), tStep(0.01), thermostat(THERMOSTAT_NONE), neighborUpdateListener(nullptr) {
  takesForces = true;
}

IntegratorMD::~IntegratorMD(){
  delete neighborUpdateListener;
}

void IntegratorMD::setNbrCheckInterval(int i) {
  if (!neighborUpdateListener && i>0) {
    neighborUpdateListener = new NeighborUpdateListener(static_cast<PotentialMasterList&>(potentialMaster), i);
    addListener(neighborUpdateListener);
    return;
  }
  if (i==0) {
    removeListener(neighborUpdateListener);
    delete neighborUpdateListener;
    return;
  }
  neighborUpdateListener->setInterval(i);
}

void IntegratorMD::allComputeFinished(double uTotNew, double virialTotNew, double **f, double* virialTensor) {
  energy = uTotNew;
  forces = f;
}

void IntegratorMD::setTimeStep(double t) {
  tStep = t;
}

void IntegratorMD::randomizeVelocity(int iAtom) {
  int iType = box.getAtomType(iAtom);
  double m = atomInfo.getMass(iType);
  double sqrtTM = sqrt(temperature/m);
  double* vi = box.getAtomVelocity(iAtom);
  for (int j=0; j<3; j++) {
    kineticEnergy -= 0.5*m*vi[j]*vi[j];
    vi[j] = random.nextGaussian()*sqrtTM;
    kineticEnergy += 0.5*m*vi[j]*vi[j];
  }
}

void IntegratorMD::randomizeVelocities(bool zeroMomentum) {
  double momentum[3] = {0};
  double totalMass = 0;
  int numAtoms = box.getNumAtoms();
  double sqrtTM[atomInfo.getNumTypes()], imass[atomInfo.getNumTypes()];
  for (int i=0; i<atomInfo.getNumTypes(); i++) {
    imass[i] = atomInfo.getMass(i);
    sqrtTM[i] = sqrt(temperature/imass[i]);
  }
  kineticEnergy = 0;
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] = random.nextGaussian()*sqrtTM[iType];
      if (zeroMomentum) {
        momentum[j] += vi[j]*imass[iType];
        totalMass += imass[iType];
      }
      else {
        kineticEnergy += 0.5*imass[iType]*vi[j]*vi[j];
      }
    }
  }
  if (zeroMomentum) {
    kineticEnergy = 0;
    for (int j=0; j<3; j++ ){
      momentum[j] /= totalMass;
    }
    for (int iAtom=0; iAtom<numAtoms; iAtom++) {
      int iType = box.getAtomType(iAtom);
      double* vi = box.getAtomVelocity(iAtom);
      for (int j=0; j<3; j++ ){
        vi[j] -= momentum[j];
        kineticEnergy += 0.5*imass[iType]*vi[j]*vi[j];
      }
    }
  }
}

void IntegratorMD::reset() {
  potentialMaster.init();
  Integrator::reset();

  kineticEnergy = 0;
  int numAtoms = box.getNumAtoms();
  double imass[atomInfo.getNumTypes()];
  for (int i=0; i<atomInfo.getNumTypes(); i++) {
    imass[i] = atomInfo.getMass(i);
  }
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      kineticEnergy += 0.5*imass[iType]*vi[j]*vi[j];
    }
  }
}

double IntegratorMD::getKineticEnergy() {
  return kineticEnergy;
}

void IntegratorMD::addPotentialCallback(PotentialCallback* callback, int interval) {
  struct PotentialCallbackInfo pci;
  pci.pcb = callback;
  pci.interval = pci.countdown = interval;
  allPotentialCallbacks.push_back(pci);
}

void IntegratorMD::computeForces() {
  selfPotentialCallbackVec.clear();
  selfPotentialCallbackVec.push_back(this);
  for (vector<struct PotentialCallbackInfo>::iterator it = allPotentialCallbacks.begin(); it!=allPotentialCallbacks.end(); it++) {
    (*it).countdown--;
    if ((*it).countdown == 0) {
      (*it).countdown = (*it).interval;
      (*it).pcb->reset();
      selfPotentialCallbackVec.push_back((*it).pcb);
    }
  }
  potentialMaster.computeAll(selfPotentialCallbackVec);
}

double** IntegratorMD::getForces() {
  return forces;
}

double IntegratorMD::getTimeStep() {
  return tStep;
}
