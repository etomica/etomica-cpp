#include <algorithm>
#include <stdio.h>
#include "potential-master.h"
#include "alloc2d.h"

PotentialCallback::PotentialCallback() : callPair(false), callFinished(false), takesForces(false) {}

PotentialMaster::PotentialMaster(SpeciesList& sl, Box& b) : speciesList(sl), box(b), force(nullptr), pureAtoms(false) {
  numAtomTypes = sl.getAtomInfo().getNumTypes();

  pairPotentials = (Potential***)malloc2D(numAtomTypes, numAtomTypes, sizeof(Potential*));
  pairCutoffs = (double**)malloc2D(numAtomTypes, numAtomTypes, sizeof(double));
  for (int i=0; i<numAtomTypes; i++) {
    for (int j=0; j<numAtomTypes; j++) {
      pairPotentials[i][j] = nullptr;
      pairCutoffs[i][j] = 0;
    }
  }
  uAtom.resize(b.getNumAtoms());
  bondedPairs = new vector<vector<int*>>[sl.size()];
  bondedPotentials = new vector<Potential*>[sl.size()];
  bondedAtoms = new vector<int>*[sl.size()];
}

void PotentialMaster::setPairPotential(int iType, int jType, Potential* p) {
  for (int i=0; i<numAtomTypes; i++) {
    for (int j=0; j<numAtomTypes; j++) {
      pairPotentials[i][j] = p;
      double rc = p->getCutoff();
      pairCutoffs[i][j] = rc*rc;
    }
  }
}

void PotentialMaster::setBondPotential(int iSpecies, vector<int*> &bp, Potential *p) {
  pureAtoms = false;
  bondedPairs[iSpecies].push_back(bp);
  bondedPotentials[iSpecies].push_back(p);
  Species* s = speciesList.get(iSpecies);
  // bondedAtoms keeps track of all atoms bonded to a given atom, so that they can be
  // excluded as a pair for standard (LJ) interactions
  bondedAtoms[iSpecies] = new vector<int>[s->getNumAtoms()];
  vector<int> *&myBA = bondedAtoms[iSpecies];
  for (int i=0; i<(int)bp.size(); i++) {
    int a0 = bp[i][0];
    int a1 = bp[i][1];
    if (!binary_search(myBA[a0].begin(), myBA[a0].end(), a1)) myBA[a0].push_back(a1);
    if (!binary_search(myBA[a1].begin(), myBA[a1].end(), a0)) myBA[a1].push_back(a0);
  }
  for (int i=0; i<s->getNumAtoms(); i++) {
    sort(myBA[i].begin(), myBA[i].end());
  }
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
  double u, du, d2u;
  double dr[3];
  for (int i=0; i<numAtoms; i++) {
    uAtom[i] = 0;
    if (doForces) for (int k=0; k<3; k++) force[i][k] = 0;
  }
  for (int i=0; i<numAtoms; i++) {
    double *ri = box.getAtomPosition(i);
    int iType = box.getAtomType(i);
    for (int j=i+1; j<numAtoms; j++) {
      int jType = box.getAtomType(j);
      double *rj = box.getAtomPosition(j);
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      box.nearestImage(dr);
      double r2 = 0;
      for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
      double rc2 = pairCutoffs[iType][jType];
      if (r2 > rc2) continue;
      pairPotentials[iType][jType]->u012(r2, u, du, d2u);
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
  uAtomsChanged.resize(0);
}

void PotentialMaster::processAtomU(int coeff) {
  int numAtomsChanged = uAtomsChanged.size();
  for (int i=0; i<numAtomsChanged; i++) {
    int iAtom = uAtomsChanged[i];
    uAtom[iAtom] += coeff*duAtom[i];
  }
  uAtomsChanged.resize(0);
}

void PotentialMaster::computeOne(int iAtom, double *ri, double &u1, bool isTrial) {
  int numAtoms = box.getNumAtoms();
  u1 = 0;
  double dr[3];
  uAtomsChanged.resize(1);
  duAtom.resize(1);
  uAtomsChanged[0] = iAtom;
  duAtom[0] = 0;
  int iType = box.getAtomType(iAtom);
  double* iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];
  for (int j=0; j<numAtoms; j++) {
    if (j==iAtom) continue;
    int jType = box.getAtomType(j);
    double *rj = box.getAtomPosition(j);
    for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
    box.nearestImage(dr);
    double r2 = 0;
    for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
    if (r2 > iCutoffs[jType]) continue;
    uAtomsChanged.push_back(j);
    double uij = iPotentials[jType]->u(r2);
    duAtom[0] += 0.5*uij;
    duAtom.push_back(0.5*uij);
    u1 += uij;
  }
}

void PotentialMaster::newAtom() {
  int n = box.getNumAtoms();
  uAtom.resize(n);
  uAtom[n-1] = 0;
  duAtom.resize(n);
}

void PotentialMaster::removeAtom(int iAtom) {
  int n = box.getNumAtoms();
  uAtom[iAtom] = uAtom[n-1];
  uAtom.resize(n);
  duAtom.resize(n);
}
