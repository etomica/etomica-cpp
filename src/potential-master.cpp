#include <algorithm>
#include <stdio.h>
#include "potential-master.h"
#include "alloc2d.h"

PotentialCallback::PotentialCallback() : callPair(false), callFinished(false), takesForces(false) {}

PotentialMaster::PotentialMaster(const SpeciesList& sl, Box& b) : speciesList(sl), box(b), duAtomSingle(false), duAtomMulti(false), force(nullptr), numAtomTypes(sl.getNumAtomTypes()), pureAtoms(sl.isPurelyAtomic()), rigidMolecules(true), doTruncationCorrection(true), doSingleTruncationCorrection(false) {

  pairPotentials = (Potential***)malloc2D(numAtomTypes, numAtomTypes, sizeof(Potential*));
  pairCutoffs = (double**)malloc2D(numAtomTypes, numAtomTypes, sizeof(double));
  for (int i=0; i<numAtomTypes; i++) {
    for (int j=0; j<numAtomTypes; j++) {
      pairPotentials[i][j] = nullptr;
      pairCutoffs[i][j] = 0;
    }
  }
  uAtom.resize(b.getNumAtoms());
  bondedPairs = new vector<vector<int*> >[sl.size()];
  bondedPotentials = new vector<Potential*>[sl.size()];
  bondedAtoms = new vector<int>*[sl.size()];

  numAtomsByType = new int[numAtomTypes];
  for (int i=0; i<numAtomTypes; i++) {
    numAtomsByType[i] = 0;
  }
  for (int i=0; i<sl.size(); i++) {
    Species *s = sl.get(i);
    int iNumMolecules = box.getNumMolecules(i);
    int iNumAtoms = s->getNumAtoms();
    int* atomTypes1 = s->getAtomTypes();
    for (int j=0; j<iNumAtoms; j++) {
      numAtomsByType[atomTypes1[j]] += iNumMolecules;
    }
  }
}

PotentialMaster::~PotentialMaster() {
  free2D((void**)pairPotentials);
  free2D((void**)pairCutoffs);
  free2D((void**)force);
  delete[] bondedPairs;
  delete[] bondedPotentials;
  delete[] bondedAtoms;
  delete[] numAtomsByType;
}

void PotentialMaster::setDoTruncationCorrection(bool doCorrection) {
  doTruncationCorrection = doCorrection;
  if (!doCorrection) doSingleTruncationCorrection = false;
}

void PotentialMaster::setDoSingleTruncationCorrection(bool doCorrection) {
  doSingleTruncationCorrection = doCorrection;
}

void PotentialMaster::setPairPotential(int iType, int jType, Potential* p) {
  pairPotentials[iType][jType] = pairPotentials[jType][iType] = p;
  double rc = p->getCutoff();
  pairCutoffs[iType][jType] = pairCutoffs[jType][iType] = rc*rc;
}

void PotentialMaster::setBondPotential(int iSpecies, vector<int*> &bp, Potential *p) {
  rigidMolecules = false;
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
    int iMolecule = i, iFirstAtom = i, iSpecies = 0;
    vector<int> *iBondedAtoms = nullptr;
    if (!pureAtoms) {
      box.getMoleculeInfoAtom(i, iMolecule, iSpecies, iFirstAtom);
      if (!rigidMolecules) {
        iBondedAtoms = &bondedAtoms[iSpecies][i-iFirstAtom];
      }
    }
    double *ri = box.getAtomPosition(i);
    int iType = box.getAtomType(i);
    for (int j=i+1; j<numAtoms; j++) {
      if (checkSkip(j, iSpecies, iMolecule, iBondedAtoms)) continue;
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
  if (!pureAtoms && !rigidMolecules) {
    computeAllBonds(doForces, uTot, virialTot);
  }
  computeAllTruncationCorrection(uTot, virialTot);
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}

void PotentialMaster::computeAllTruncationCorrection(double &uTot, double &virialTot) {
  if (!doTruncationCorrection) return;
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  for (int i=0; i<numAtomTypes; i++) {
    int iNumAtoms = numAtomsByType[i];
    for (int j=i; j<numAtomTypes; j++) {
      int jNumAtoms = numAtomsByType[j];
      Potential *p = pairPotentials[i][j];
      double u, du, d2u;
      p->u012TC(u, du, d2u);
      int numPairs = j==i ? iNumAtoms*(jNumAtoms-1)/2 : iNumAtoms*jNumAtoms;
      uTot += u * numPairs/vol;
      virialTot += du * numPairs/vol;
    }
  }
}

double PotentialMaster::computeOneTruncationCorrection(int iAtom) {
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  int iType = box.getAtomType(iAtom);
  double uTot = 0;
  for (int j=0; j<numAtomTypes; j++) {
    int jNumAtoms = numAtomsByType[j];
    Potential *p = pairPotentials[iType][j];
    double u, du, d2u;
    p->u012TC(u, du, d2u);
    uTot += u * jNumAtoms/vol;
  }
  return uTot;
}

void PotentialMaster::computeAllBonds(bool doForces, double &uTot, double &virialTot) {
  int numSpecies = speciesList.size();
  for (int iSpecies=0; iSpecies<numSpecies; iSpecies++) {
    vector<Potential*> &iBondedPotentials = bondedPotentials[iSpecies];
    if (iBondedPotentials.size() == 0) continue;
    int sna = speciesList.get(iSpecies)->getNumAtoms();
    vector<vector<int*> > iBondedPairs = bondedPairs[iSpecies];
    for (int j=0; j<(int)iBondedPotentials.size(); j++) {
      Potential* p = iBondedPotentials[j];
      vector<int*> jBondedPairs = iBondedPairs[j];
      int firstAtom = box.getFirstAtom(iSpecies, 0);
      for (int kMolecule=0; kMolecule<box.getNumMolecules(iSpecies); kMolecule++) {
        for (int l=0; l<(int)jBondedPairs.size(); l++) {
          int* lBondedPair = jBondedPairs[l];
          int iAtom = firstAtom + lBondedPair[0];
          int jAtom = firstAtom + lBondedPair[1];

          double *ri = box.getAtomPosition(iAtom);
          double *rj = box.getAtomPosition(jAtom);
          double dr[3];
          for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
          box.nearestImage(dr);
          double r2 = 0;
          for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
          double u, du, d2u;
          p->u012(r2, u, du, d2u);
          uAtom[iAtom] += 0.5*u;
          uAtom[jAtom] += 0.5*u;
          uTot += u;
          virialTot += du;
          for (vector<PotentialCallback*>::iterator it = pairCallbacks.begin(); it!=pairCallbacks.end(); it++) {
            (*it)->pairCompute(iAtom, jAtom, dr, u, du, d2u);
          }

          // f0 = dr du / r^2
          if (!doForces) continue;
          du /= r2;
          for (int m=0; m<3; m++) {
            force[iAtom][m] += dr[m]*du;
            force[jAtom][m] -= dr[m]*du;
          }
        }
        firstAtom += sna;
      }
    }
  }
}

void PotentialMaster::computeOneMoleculeBonds(const int iSpecies, const int iMolecule, double &u1) {
  vector<Potential*> &iBondedPotentials = bondedPotentials[iSpecies];
  if (iBondedPotentials.size() == 0) return;
  int sna = speciesList.get(iSpecies)->getNumAtoms();
  vector<vector<int*> > iBondedPairs = bondedPairs[iSpecies];
  for (int j=0; j<(int)iBondedPotentials.size(); j++) {
    Potential* p = iBondedPotentials[j];
    vector<int*> jBondedPairs = iBondedPairs[j];
    int firstAtom = box.getFirstAtom(iSpecies, iMolecule);
    for (int l=0; l<(int)jBondedPairs.size(); l++) {
      int* lBondedPair = jBondedPairs[l];
      int iAtom = firstAtom + lBondedPair[0];
      int jAtom = firstAtom + lBondedPair[1];

      double *ri = box.getAtomPosition(iAtom);
      double *rj = box.getAtomPosition(jAtom);
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      box.nearestImage(dr);
      double r2 = 0;
      for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
      double u = p->u(r2);
      uAtom[iAtom] += 0.5*u;
      uAtom[jAtom] += 0.5*u;
      u1 += u;
    }
    firstAtom += sna;
  }
}

double PotentialMaster::oldEnergy(int iAtom) {
  double u = 2*uAtom[iAtom];
  if (doSingleTruncationCorrection) {
    u += computeOneTruncationCorrection(iAtom);
  }
  return u;
}

double PotentialMaster::oldMoleculeEnergy(int iMolecule) {
  // only works for rigid molecules
  int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
  double u = 0;
  for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
    u += 2*uAtom[iAtom];
    if (doSingleTruncationCorrection) {
      u += computeOneTruncationCorrection(iAtom);
    }
  }
  return u;
}

void PotentialMaster::resetAtomDU() {
  int numAtomsChanged = uAtomsChanged.size();
  if (duAtomMulti) {
    for (int i=0; i<numAtomsChanged; i++) {
      int iAtom = uAtomsChanged[i];
      duAtom[iAtom] = 0;
    }
  }
  else {
    for (int i=0; i<numAtomsChanged; i++) duAtom[i] = 0;
  }
  uAtomsChanged.resize(0);
  duAtomSingle = duAtomMulti = false;
}

void PotentialMaster::processAtomU(int coeff) {
  if (duAtomSingle && duAtomMulti) {
    fprintf(stderr, "Can't simultaneously do single and multi duAtom!\n");
    abort();
  }
  int numAtomsChanged = uAtomsChanged.size();
  if (duAtomSingle) {
    for (int i=0; i<numAtomsChanged; i++) {
      int iAtom = uAtomsChanged[i];
      uAtom[iAtom] += coeff*duAtom[i];
      duAtom[i] = 0;
    }
  }
  else if (duAtomMulti) {
    for (int i=0; i<numAtomsChanged; i++) {
      int iAtom = uAtomsChanged[i];
      uAtom[iAtom] += coeff*duAtom[iAtom];
      duAtom[iAtom] = 0;
    }
  }
  uAtomsChanged.resize(0);
  duAtomSingle = duAtomMulti = false;
}

void PotentialMaster::computeOne(const int iAtom, double &u1) {
  duAtomSingle = true;
  u1 = 0;
  uAtomsChanged.resize(1);
  duAtom.resize(1);
  uAtomsChanged[0] = iAtom;
  duAtom[0] = 0;
  int iMolecule = 0, iFirstAtom = 0, iSpecies = 0;
  if (!pureAtoms && !rigidMolecules) {
    box.getMoleculeInfoAtom(iAtom, iMolecule, iSpecies, iFirstAtom);
  }
  const double *ri = box.getAtomPosition(iAtom);
  computeOneInternal(iAtom, ri, u1, iSpecies, iMolecule, iFirstAtom);
  if (doSingleTruncationCorrection) {
    u1 += computeOneTruncationCorrection(iAtom);
  }
}

void PotentialMaster::computeOneInternal(const int iAtom, const double *ri, double &u1, const int iSpecies, const int iMolecule, const int iFirstAtom) {
  vector<int> *iBondedAtoms = nullptr;
  if (!pureAtoms && !rigidMolecules) {
    iBondedAtoms = &bondedAtoms[iSpecies][iAtom-iFirstAtom];
  }
  int iType = box.getAtomType(iAtom);
  double* iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];
  int numAtoms = box.getNumAtoms();
  for (int jAtom=0; jAtom<numAtoms; jAtom++) {
    if (jAtom==iAtom) continue;
    if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
    int jType = box.getAtomType(jAtom);
    double *rj = box.getAtomPosition(jAtom);
    double dr[3];
    for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
    box.nearestImage(dr);
    double r2 = 0;
    for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
    if (r2 > iCutoffs[jType]) continue;
    double uij = iPotentials[jType]->u(r2);
    if (duAtomSingle) {
      uAtomsChanged.push_back(jAtom);
      duAtom[0] += 0.5*uij;
      duAtom.push_back(0.5*uij);
    }
    else {
      if (duAtom[jAtom] == 0) uAtomsChanged.push_back(jAtom);
      duAtom[iAtom] += 0.5*uij;
      duAtom[jAtom] += 0.5*uij;
    }
    u1 += uij;
  }
}

void PotentialMaster::computeOneMolecule(int iMolecule, double &u1) {
  duAtomMulti = true;
  int numAtoms = box.getNumAtoms();
  u1 = 0;
  uAtomsChanged.resize(0);
  if ((int)duAtom.size() < numAtoms) {
    int s = duAtom.size();
    duAtom.resize(numAtoms);
    fill(duAtom.begin()+s, duAtom.end(), 0);
  }
  int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
  for (int iAtom=firstAtom; iAtom<=lastAtom; iAtom++) {
    if (duAtom[iAtom] == 0) {
      uAtomsChanged.push_back(iAtom);
    }
    double *ri = box.getAtomPosition(iAtom);
    computeOneInternal(iAtom, ri, u1, iSpecies, iMolecule, firstAtom);
    if (doSingleTruncationCorrection) {
      u1 += computeOneTruncationCorrection(iAtom);
    }
  }
  if (!pureAtoms && !rigidMolecules) {
    computeOneMoleculeBonds(iSpecies, iMolecule, u1);
  }
}

void PotentialMaster::newMolecule(int iSpecies) {
  int iMolecule = box.getNumMolecules(iSpecies)-1;
  int firstAtom = box.getFirstAtom(iSpecies, iMolecule);
  int speciesAtoms = speciesList.get(iSpecies)->getNumAtoms();
  int lastAtom = firstAtom + speciesAtoms - 1;
  int numAtoms = box.getNumAtoms();
  uAtom.resize(numAtoms);
  // first we have to shift uAtoms for all species>iSpecies
  for (int jAtom=numAtoms-1; jAtom>lastAtom; jAtom--) {
    uAtom[jAtom] = uAtom[jAtom-speciesAtoms];
  }
  for (int jAtom=lastAtom; jAtom>=firstAtom; jAtom--) {
    uAtom[jAtom] = 0;
    numAtomsByType[box.getAtomType(jAtom)]++;
  }
}

void PotentialMaster::removeMolecule(int iSpecies, int iMolecule) {
  int firstAtom = box.getFirstAtom(iSpecies, iMolecule);
  int speciesAtoms = speciesList.get(iSpecies)->getNumAtoms();
  int jMolecule = box.getNumMolecules(iSpecies)-1;
  int jFirstAtom = box.getFirstAtom(iSpecies, jMolecule);
  for (int i=0; i<speciesAtoms; i++) {
    uAtom[firstAtom+i] = uAtom[jFirstAtom+i];
    numAtomsByType[box.getAtomType(firstAtom+i)]--;
  }
  int numAtoms = box.getNumAtoms();
  // now shift uAtoms for all species>iSpecies
  for (int jAtom=jFirstAtom+speciesAtoms; jAtom<numAtoms; jAtom++) {
    uAtom[jAtom-speciesAtoms] = uAtom[jAtom];
  }
  uAtom.resize(numAtoms-speciesAtoms);
}

double PotentialMaster::uTotalFromAtoms() {
  double uTot = 0;
  int numAtoms = box.getNumAtoms();
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    uTot += uAtom[iAtom];
  }
  double virialTot;
  computeAllTruncationCorrection(uTot, virialTot);
  return uTot;
}
