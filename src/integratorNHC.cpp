#include "integrator.h"

/**
 * Nose-Hoover chains molecular dynamics integrator
 * adopted from Allen & Tildesley's example code, md_nvt_lj
 */

IntegratorNHC::IntegratorNHC(AtomInfo& ai, PotentialMaster& p, Random& r, Box& b, int nc, double tau) : IntegratorMD(ai, p, r, b), numChains(nc) {
  q = new double[numChains];
  eta = new double[numChains];
  etaP = new double[numChains];
  double g = 3*(box.getNumAtoms()-1);
  q[0] = g*temperature*tau*tau;
  for (int i=1; i<numChains; i++) q[i] = temperature * tau*tau;
}

IntegratorNHC::~IntegratorNHC(){
  delete[] q;
  delete[] eta;
  delete[] etaP;
}

void IntegratorNHC::setTemperature(double newTemperature) {
  IntegratorMD::setTemperature(newTemperature);
}

void IntegratorNHC::doStep() {
  stepCount++;
  propagatorU4(tStep/4, -1);
  propagatorU3(tStep/2);
  propagatorU4(tStep/4, +1);
  propagatorU2(tStep/2);
  propagatorU1(tStep);

  computeForces();

  propagatorU2(tStep/2);
  propagatorU4(tStep/4, -1);
  propagatorU3(tStep/2);
  propagatorU4(tStep/4, +1);

  for (vector<IntegratorListener*>::iterator it = listenersStepFinished.begin(); it!=listenersStepFinished.end(); it++) {
    (*it)->stepFinished();
  }
  if (nbrCheckInterval>0) nbrCheckCountdown--;
#ifdef DEBUG
  printf("%ld step %e %e %e\n", stepCount, energy/n, kineticEnergy/n, (energy + kineticEnergy)/n);
#endif
}

void IntegratorNHC::propagatorU1(double dt) {
  int n = box.getNumAtoms();
  for (int iAtom=0; iAtom<n; iAtom++) {
    double* ri = box.getAtomPosition(iAtom);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      ri[j] += dt*vi[j];
    }
  }
}

void IntegratorNHC::propagatorU2(double dt) {
  int n = box.getNumAtoms();
  kineticEnergy = 0;
  for (int iAtom=0; iAtom<n; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double mass = atomInfo.getMass(iType);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] += dt*forces[iAtom][j]/mass;
      kineticEnergy += 0.5*mass*vi[j]*vi[j];
    }
  }
}

void IntegratorNHC::propagatorU3(double dt) {
  int n = box.getNumAtoms();
  double fac = exp(-dt*etaP[0]/q[0]);
  //printf("fac %e  %e %e\n", fac, etaP[0], q[0]);
  kineticEnergy = 0;
  for (int iAtom=0; iAtom<n; iAtom++) {
    int iType = box.getAtomType(iAtom);
    double mass = atomInfo.getMass(iType);
    double* vi = box.getAtomVelocity(iAtom);
    for (int j=0; j<3; j++) {
      vi[j] *= fac;
      kineticEnergy += 0.5*mass*vi[j]*vi[j];
    }
  }
  for (int i=0; i<numChains; i++) {
    eta[i] += dt*etaP[i]/q[i];
  }
}

void IntegratorNHC::propagatorU4(double dt, int direction) {
  double g = 3*(box.getNumAtoms()-1);
  //printf("%f %f\n", temperature, kineticEnergy*2/g);
  for (int i=0; i<numChains; i++) {
    int ic = direction == 1 ? i : (numChains-i-1);
    double gi = ic==0 ? (2*kineticEnergy - g*temperature) : ((etaP[ic-1]*etaP[ic-1]/q[ic-1]) - temperature);
    if (ic == numChains-1) {
      etaP[ic] += dt*gi;
      //printf("etaP[0] %e\n", etaP[0]);
      continue;
    }
    double x = dt * etaP[ic+1]/q[ic+1];
    double c;
    if (x < 0.001) {
      // for small values, use Taylor series
      c = 1 + x*(-0.5 + x*(1.0/6.0 - x/24));
    }
    else {
      c = (1-exp(-x))/x;
    }
    etaP[ic] = etaP[ic] * exp(-x) + dt*gi*c;
  }
}

void IntegratorNHC::reset() {
  IntegratorMD::reset();
  randomizeVelocities(true);
}
