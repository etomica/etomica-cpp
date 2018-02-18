#include "integrator.h"

IntegratorMD::IntegratorMD(PotentialMaster& p, Random& r, Box& b) : Integrator(p), random(r), box(b), forces(NULL), tStep(0.01), thermostat(THERMOSTAT_NONE), nbrCheckInterval(-1), nbrCheckCountdown(-1) {
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

void IntegratorMD::doStep() {
  if (nbrCheckCountdown==0) {
    static_cast<PotentialMasterList&>(potentialMaster).checkUpdateNbrs();
    nbrCheckCountdown = nbrCheckInterval;
  }
  stepCount++;
  int n = box.getNumAtoms();
  for (int iAtom=0; iAtom<n; iAtom++) {
    double* ri = box.getAtomPosition(iAtom);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] += 0.5*tStep*forces[iAtom][j];
      ri[j] += tStep*vi[j];
    }
  }

  potentialMaster.computeAll(selfPotentialCallbackVec);

  for (int iAtom=0; iAtom<n; iAtom++) {
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] += 0.5*tStep*forces[iAtom][j];
    }
  }
  for (vector<IntegratorListener*>::iterator it = listeners.begin(); it!=listeners.end(); it++) {
    (*it)->stepFinished();
  }
  if (nbrCheckInterval>0) nbrCheckCountdown--;
}

void IntegratorMD::randomizeVelocities(bool zeroMomentum) {
  double imass = 1;
  double sqrtTM = sqrt(temperature/imass);
  double momentum[3];
  double totalMass = 0;
  int numAtoms = box.getNumAtoms();
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] = random.nextGaussian()*sqrtTM;
      if (zeroMomentum) {
        momentum[j] += vi[j]*imass;
        totalMass += imass;
      }
    }
  }
  if (zeroMomentum) {
    for (int j=0; j<3; j++ ){
      momentum[j] /= totalMass;
    }
    for (int iAtom=0; iAtom<numAtoms; iAtom++) {
      double* vi = box.getAtomVelocity(iAtom);
      for (int j=0; j<3; j++ ){
        vi[j] -= momentum[j];
      }
    }
  }
}

void IntegratorMD::reset() {
  Integrator::reset();
  randomizeVelocities(false);
}
