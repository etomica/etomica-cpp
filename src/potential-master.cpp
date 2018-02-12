#include <stdio.h>
#include "potential-master.h"
#include "alloc2d.h"

PotentialCallback::PotentialCallback() : callPair(false), callFinished(false), takesForces(false) {}

PotentialMaster::PotentialMaster(Potential& p2, Box& b) : potential(p2), box(b) {
  uAtom.resize(b.getNumAtoms());
  duAtom.resize(b.getNumAtoms());
  uAtomsChanged.resize(b.getNumAtoms());
  numAtomsChanged = 0;
  force = nullptr;
}

Box& PotentialMaster::getBox() {
  return box;
}

void PotentialMaster::computeAll(vector<PotentialCallback*> &callbacks) {
  pairCallbacks.resize(0);
  bool doForces = false;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callPair) pairCallbacks.push_back(*it);
    if ((*it)->takesForces) doForces = true;
  }
  if (doForces && !force) {
    force = (double**)malloc2D(box.getNumAtoms(), 3, sizeof(double));
  }

  int numAtoms = box.getNumAtoms();
  double uTot = 0, virialTot = 0;
  double rc = potential.getCutoff();
  double rc2 = rc*rc;
  double u, du, d2u;
  double dr[3];
  for (int i=0; i<numAtoms; i++) {
    uAtom[i] = 0;
    if (doForces) for (int k=0; k<3; k++) force[i][k] = 0;
  }
  for (int i=0; i<numAtoms; i++) {
    double *ri = box.getAtomPosition(i);
    for (int j=i+1; j<numAtoms; j++) {
      double *rj = box.getAtomPosition(j);
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      box.nearestImage(dr);
      double r2 = 0;
      for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
      if (r2 > rc2) continue;
      potential.u012(r2, u, du, d2u);
      uAtom[i] += 0.5*u;
      uAtom[j] += 0.5*u;
      uTot += u;
      virialTot += du;
      for (vector<PotentialCallback*>::iterator it = pairCallbacks.begin(); it!=pairCallbacks.end(); it++) {
        (*it)->pairCompute(i, j, dr, u, du, d2u);
      }

      // f0 = dr du / r^2
      if (!doForces) continue;
      du /= r2;
      for (int k=0; k<3; k++) {
        force[i][k] += dr[k]*du;
        force[j][k] -= dr[k]*du;
      }
    }
  }
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}

double PotentialMaster::oldEnergy(int iAtom) {
  return 2*uAtom[iAtom];
}

void PotentialMaster::resetAtomDU() {
  numAtomsChanged = 0;
}

void PotentialMaster::processAtomU(int coeff) {
  for (int i=0; i<numAtomsChanged; i++) {
    int iAtom = uAtomsChanged[i];
    uAtom[iAtom] += coeff*duAtom[i];
  }
}

void PotentialMaster::computeOne(int iAtom, double *ri, double &u1, bool isTrial) {
  int numAtoms = box.getNumAtoms();
  u1 = 0;
  double rc = potential.getCutoff();
  double rc2 = rc*rc;
  double dr[3];
  numAtomsChanged = 1;
  uAtomsChanged[0] = iAtom;
  duAtom[0] = 0;
  for (int j=0; j<numAtoms; j++) {
    if (j==iAtom) continue;
    double *rj = box.getAtomPosition(j);
    for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
    box.nearestImage(dr);
    double r2 = 0;
    for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
    if (r2 > rc2) continue;
    uAtomsChanged[numAtomsChanged] = j;
    double uij = potential.u(r2);
    duAtom[0] += 0.5*uij;
    duAtom[numAtomsChanged] = 0.5*uij;
    numAtomsChanged++;
    u1 += uij;
  }
}

void PotentialMaster::newAtom() {
  int n = box.getNumAtoms();
  uAtom.resize(n);
  uAtom[n-1] = 0;
  duAtom.resize(n);
  uAtomsChanged.resize(n);
}

void PotentialMaster::removeAtom(int iAtom) {
  int n = box.getNumAtoms();
  uAtom[iAtom] = uAtom[n-1];
  uAtom.resize(n);
  duAtom.resize(n);
  uAtomsChanged.resize(n);
}
