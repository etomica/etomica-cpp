#include "integrator.h"

IntegratorNVE::IntegratorNVE(AtomInfo& ai, PotentialMaster& p, Random& r, Box& b) : IntegratorMD(ai, p, r, b) {
}

IntegratorNVE::~IntegratorNVE(){}

void IntegratorNVE::doStep() {
  if (nbrCheckCountdown==0) {
    static_cast<PotentialMasterList&>(potentialMaster).checkUpdateNbrs();
    nbrCheckCountdown = nbrCheckInterval;
  }
  stepCount++;
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
  if (nbrCheckInterval>0) nbrCheckCountdown--;
#ifdef DEBUG
  printf("%ld step %e %e %e\n", stepCount, energy/n, kineticEnergy/n, (energy + kineticEnergy)/n);
#endif
}

void IntegratorNVE::reset() {
  IntegratorMD::reset();
  randomizeVelocities(false);
}
