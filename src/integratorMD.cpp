#include "integrator.h"
#include "potential-master.h"

IntegratorMD::IntegratorMD(AtomInfo& ai, PotentialMaster& p, Random& r, Box& b) : Integrator(p), atomInfo(ai), random(r), box(b), forces(NULL), tStep(0.01), thermostat(THERMOSTAT_NONE), nbrCheckInterval(-1), nbrCheckCountdown(-1) {
  takesForces = true;
}

IntegratorMD::~IntegratorMD(){}

void IntegratorMD::setNbrCheckInterval(int i) {
  nbrCheckInterval = nbrCheckCountdown = i;
}

void IntegratorMD::allComputeFinished(double uTotNew, double virialTotNew, double **f) {
  energy = uTotNew;
  forces = f;
}

void IntegratorMD::setTimeStep(double t) {
  tStep = t;
}

void IntegratorMD::randomizeVelocities(bool zeroMomentum) {
  double momentum[3];
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
  if (nbrCheckInterval>0) {
    static_cast<PotentialMasterList&>(potentialMaster).reset();
    nbrCheckCountdown = nbrCheckInterval;
  }
  Integrator::reset();
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
