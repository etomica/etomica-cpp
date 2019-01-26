/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <algorithm>
#include <stdio.h>
#include "potential-master.h"
#include "potential-angle.h"
#include "ewald.h"
#include "alloc2d.h"
#include "util.h"

PotentialCallback::PotentialCallback() : callPair(false), callFinished(false), takesForces(false), takesPhi(false), takesDFDV(false), takesVirialTensor(false) {}

PotentialMaster::PotentialMaster(const SpeciesList& sl, Box& b, bool doEmbed) : speciesList(sl), box(b), duAtomSingle(false), duAtomMulti(false), force(nullptr), numForceAtoms(0), numRhoSumAtoms(0), rhoSum(nullptr), idf(nullptr), numAtomTypes(sl.getNumAtomTypes()), pureAtoms(sl.isPurelyAtomic()), rigidMolecules(true), doTruncationCorrection(true), doSingleTruncationCorrection(false), embeddingPotentials(doEmbed), charges(nullptr), sFacAtom(nullptr), doEwald(false), ewald(nullptr) {

  if (embeddingPotentials && !pureAtoms) {
    fprintf(stderr, "Embedding potentials require a purely atomic system");
    abort();
  }
  pairPotentials = (Potential***)malloc2D(numAtomTypes, numAtomTypes, sizeof(Potential*));
  pairCutoffs = (double**)malloc2D(numAtomTypes, numAtomTypes, sizeof(double));
  if (embeddingPotentials) {
    rhoPotentials = (Potential**)malloc(numAtomTypes*sizeof(Potential*));
    embedF = (EmbedF**)malloc(numAtomTypes*sizeof(EmbedF*));
    rhoCutoffs = (double*)malloc(numAtomTypes*sizeof(double));
    setDoTruncationCorrection(false);
    rhoSum = (double*)malloc(b.getNumAtoms()*sizeof(double));
    drhoSum.resize(b.getNumAtoms());
    fill(drhoSum.begin(), drhoSum.end(), 0);
    numRhoSumAtoms = b.getNumAtoms();
  }
  else {
    rhoPotentials = nullptr;
    embedF = nullptr;
    rhoCutoffs = nullptr;
  }
  for (int i=0; i<numAtomTypes; i++) {
    for (int j=0; j<numAtomTypes; j++) {
      pairPotentials[i][j] = nullptr;
      pairCutoffs[i][j] = 0;
    }
    if (embeddingPotentials) {
      rhoPotentials[i] = nullptr;
      rhoCutoffs[i] = 0;
      embedF[i] = nullptr;
    }
  }
  bondedPairs = new vector<vector<int*> >[sl.size()];
  bondedPotentials = new vector<Potential*>[sl.size()];
  bondedAtoms = new vector<int>*[sl.size()];
  bondAngleTriplets = new vector<vector<int*> >[sl.size()];
  bondAnglePotentials = new vector<PotentialAngle*>[sl.size()];

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
  free(rhoPotentials);
  free(embedF);
  free2D((void**)pairCutoffs);
  free(rhoCutoffs);
  free(rhoSum);
  free2D((void**)force);
  free(idf);
  delete[] bondedPairs;
  delete[] bondedPotentials;
  delete[] bondedAtoms;
  delete[] bondAngleTriplets;
  delete[] bondAnglePotentials;
  delete[] numAtomsByType;
}

void PotentialMaster::init() {
  int na = box.getNumAtoms();
  int s = uAtom.size();
  if (s < na || s > 2*na) {
    uAtom.resize(na);
  }
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
  if (!p) {
    pairCutoffs[iType][jType] = 0;
    return;
  }
  double rc = p->getCutoff();
  pairCutoffs[iType][jType] = pairCutoffs[jType][iType] = rc*rc;
}

void PotentialMaster::setRhoPotential(int jType, Potential* p) {
  if (!embeddingPotentials) {
    fprintf(stderr, "Potential master not configured for embedding potentials!\n");
    abort();
  }
  rhoPotentials[jType] = p;
  double rc = p->getCutoff();
  rhoCutoffs[jType] = rc*rc;
}

void PotentialMaster::setEmbedF(int iType, EmbedF* ef) {
  if (!embeddingPotentials) {
    fprintf(stderr, "Potential master not configured for embedding potentials!\n");
    abort();
  }
  embedF[iType] = ef;
}

void PotentialMaster::setEwald(EwaldBase* e) {
  doEwald = e != nullptr;
  ewald = e;
  if (doEwald) ewald->setNumAtomsByType(numAtomsByType);
}

void PotentialMaster::setBondPotential(int iSpecies, vector<int*> &bp, Potential *p) {
  if (pureAtoms) {
    fprintf(stderr, "Potential master was configured for purely atomic interactions\n");
    abort();
  }
  rigidMolecules = false;
  if (ewald) {
    ewald->setRigidMolecules(false);
  }
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

void PotentialMaster::setBondAnglePotential(int iSpecies, vector<int*> &bt, PotentialAngle *p) {
  if (pureAtoms) {
    fprintf(stderr, "Potential master was configured for purely atomic interactions\n");
    abort();
  }
  rigidMolecules = false;
  bondAngleTriplets[iSpecies].push_back(bt);
  bondAnglePotentials[iSpecies].push_back(p);
  Species* s = speciesList.get(iSpecies);
  // bondedAtoms keeps track of all atoms bonded to a given atom, so that they can be
  // excluded as a pair for standard (LJ) interactions
  bondedAtoms[iSpecies] = new vector<int>[s->getNumAtoms()];
  vector<int> *&myBA = bondedAtoms[iSpecies];
  for (int i=0; i<(int)bt.size(); i++) {
    for (int j=0; j<3; j++) {
      int aj = bt[i][j];
      for (int k=j+1; k<3; k++) {
        int ak = bt[i][k];
        if (!binary_search(myBA[aj].begin(), myBA[aj].end(), ak)) myBA[aj].push_back(ak);
        if (!binary_search(myBA[ak].begin(), myBA[ak].end(), aj)) myBA[ak].push_back(aj);
      }
    }
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
  bool doForces = false, doPhi = false, doDFDV = false, doVirialTensor = false;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if (!embeddingPotentials && (*it)->callPair) pairCallbacks.push_back(*it);
    if ((*it)->takesForces) doForces = true;
    if ((*it)->takesPhi) doPhi = true;
    if ((*it)->takesDFDV) doDFDV = true;
    if ((*it)->takesVirialTensor) doVirialTensor = true;
  }
  int numAtoms = box.getNumAtoms();
  if (doForces && numAtoms > numForceAtoms) {
    force = (double**)realloc2D((void**)force, numAtoms, 3, sizeof(double));
    if (embeddingPotentials) {
      idf = (double*)realloc(idf, numAtoms*sizeof(double));
    }
    if (doEwald) {
      sFacAtom = (complex<double>*)realloc((void**)sFacAtom, box.getNumAtoms()*sizeof(double));
    }
    numForceAtoms = numAtoms;
  }
  if (embeddingPotentials && numAtoms > numRhoSumAtoms) {
    drhoSum.resize(numAtoms);
    numRhoSumAtoms = numAtoms;
    rhoSum = (double*)realloc(rhoSum, numAtoms*sizeof(double));
  }

  double uTot = 0, virialTot = 0;
  double virialTensor[6] = {0,0,0,0,0,0};
  double dr[3];
  double zero[3];
  zero[0] = zero[1] = zero[2] = 0;
  for (int i=0; i<numAtoms; i++) {
    uAtom[i] = 0;
    if (embeddingPotentials) {
      rhoSum[i] = 0;
      drhoSum[i] = 0;
    }
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
    Potential* iRhoPotential = embeddingPotentials ? rhoPotentials[iType] : nullptr;
    double iRhoCutoff = embeddingPotentials ? rhoCutoffs[iType] : 0;
    for (int j=0; j<i; j++) {
      if (checkSkip(j, iSpecies, iMolecule, iBondedAtoms)) continue;
      int jType = box.getAtomType(j);
      Potential* pij = pairPotentials[iType][jType];
      if (!pij) continue;
      double *rj = box.getAtomPosition(j);
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      box.nearestImage(dr);
      handleComputeAll(i, j, zero, dr, zero, pij, uAtom[i], uAtom[j], doForces?force[i]:nullptr, doForces?force[j]:nullptr, uTot, virialTot, virialTensor, pairCutoffs[iType][jType], iRhoPotential, iRhoCutoff, iType, jType, doForces, doVirialTensor, false);
    }
  }
  if (embeddingPotentials) {
    // we need another pass to include embedding contributions
    int rdrhoIdx = 0;
    for (int iAtom=0; iAtom<numAtoms; iAtom++) {
      int iType = box.getAtomType(iAtom);
      double f, df, d2f;
      embedF[iType]->f012(rhoSum[iAtom], f, df, d2f);
      uTot += f;
      if (doForces) {
        double *ri = box.getAtomPosition(iAtom);
        idf[iAtom] = df;
        double iRhoCutoff = rhoCutoffs[iType];
        Potential* iRhoPotential = rhoPotentials[iType];
        for (int jAtom=0; jAtom<iAtom; jAtom++) {
          int jType = box.getAtomType(jAtom);
          double *rj = box.getAtomPosition(jAtom);
          for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
          box.nearestImage(dr);

          handleComputeAllEmbed(iAtom, jAtom, iType, jType, zero, dr, zero, df, virialTot, iRhoPotential, iRhoCutoff, rdrhoIdx);
        }
      }
    }
    if (doForces) {
      rdrho.clear();
    }
  }
  if (doEwald) {
    ewald->computeAllFourier(doForces, doPhi, doDFDV, doVirialTensor, uTot, virialTot, force, virialTensor, &pairCallbacks);
  }
  if (doForces && !pureAtoms) {
    virialTot += computeVirialIntramolecular(doVirialTensor, virialTensor);
  }
  if (!pureAtoms && !rigidMolecules) {
    computeAllBonds(doForces, uTot);
  }
  computeAllTruncationCorrection(uTot, virialTot);
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force, virialTensor);
  }
}

double PotentialMaster::computeVirialIntramolecular(const bool doVirialTensor, double virialTensor[6]) {
  double virialTot = 0;
  for (int iMolecule=0; iMolecule<box.getTotalNumMolecules(); iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    if (iLastAtom==iFirstAtom) continue;
    double* center = speciesList.get(iSpecies)->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      double* ri = box.getAtomPosition(iAtom);
      double dr[3];
      for (int k=0; k<3; k++) dr[k] = ri[k]-center[k];
      box.nearestImage(dr);
      for (int k=0; k<3; k++) virialTot += force[iAtom][k] * dr[k];
      if (doVirialTensor) {
        int l = 0;
        for (int k=0; k<3; k++) {
          for (int m=k; m<3; m++) {
            virialTensor[l] += force[iAtom][k] * dr[m];
            l++;
          }
        }
      }
    }
  }
  return virialTot;
}

void PotentialMaster::computeAllTruncationCorrection(double &uTot, double &virialTot) {
  if (!doTruncationCorrection) return;
  const double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];
  for (int i=0; i<numAtomTypes; i++) {
    int iNumAtoms = numAtomsByType[i];
    for (int j=i; j<numAtomTypes; j++) {
      Potential *p = pairPotentials[i][j];
      if (p==nullptr) continue;
      int jNumAtoms = numAtomsByType[j];
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

void PotentialMaster::handleOneBondPair(bool doForces, double &uTot, int iAtom, int jAtom, Potential* p) {
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
  for (vector<PotentialCallback*>::iterator it = pairCallbacks.begin(); it!=pairCallbacks.end(); it++) {
    (*it)->pairCompute(iAtom, jAtom, dr, u, du, d2u);
  }
  if (!doForces) return;

  // f0 = dr du / r^2
  du /= r2;
  for (int m=0; m<3; m++) {
    force[iAtom][m] += dr[m]*du;
    force[jAtom][m] -= dr[m]*du;
  }
}

void PotentialMaster::handleOneBondAngleTriplet(bool doForces, double &uTot, int iAtom, int jAtom, int kAtom, PotentialAngle* p) {
  double *ri = box.getAtomPosition(iAtom);
  double *rj = box.getAtomPosition(jAtom);
  double *rk = box.getAtomPosition(kAtom);
  double drji[3], drjk[3];
  for (int l=0; l<3; l++) {
    drji[l] = ri[l]-rj[l];
    drjk[l] = rk[l]-rj[l];
  }
  box.nearestImage(drji);
  box.nearestImage(drjk);
  double rji2 = 0, rjk2 = 0;
  for (int l=0; l<3; l++) {
    rji2 += drji[l]*drji[l];
    rjk2 += drjk[l]*drjk[l];
  }
  double jidotjk = 0;
  for (int l=0; l<3; l++) {
    jidotjk += drji[l]*drjk[l];
  }
  double rji1rjk1 = sqrt(rji2*rjk2);
  double costheta = jidotjk/rji1rjk1;
  if (costheta<-1) costheta=-1;
  else if (costheta>1) costheta=1;
  double u, du, d2u;
  p->u012(costheta, u, du, d2u);
  uAtom[iAtom] += 0.3333333333333333333*u;
  uAtom[jAtom] += 0.3333333333333333333*u;
  uAtom[kAtom] += 0.3333333333333333333*u;
  uTot += u;
  for (vector<PotentialCallback*>::iterator it = pairCallbacks.begin(); it!=pairCallbacks.end(); it++) {
    //(*it)->Compute(iAtom, jAtom, dr, u, du, d2u);
  }
  if (!doForces || fabs(costheta) > 0.9999999999) return;

  // each end atom has a force perpendicular to its bond and
  // inversely proportional to its bond length
  // central atom has a force to conteract the forces on the end atoms
  double rji1 = sqrt(rji2);
  double rjk1 = rji1rjk1 / rji1;
  double fi[3], fk[3];
  double fi2 = 0, fk2 = 0;
  for (int l=0; l<3; l++) {
    fi[l] = drjk[l] - jidotjk/rji1;
    fi2 += fi[l]*fi[l];
    fk[l] = drji[l] - jidotjk/rjk1;
    fk2 += fk[l]*fk[l];
  }
  double fi1 = sqrt(fi2), fk1 = sqrt(fk2);
  for  (int l=0; l<3; l++) {
    fi[l] *= du/(fi1*rji1);
    force[iAtom][l] += fi[l];
    fk[l] *= du/(fk1*rjk1);
    force[kAtom][l] += fk[l];
    force[jAtom][l] = -(fi[l] + fk[l]);
  }
}

void PotentialMaster::computeAllBonds(bool doForces, double &uTot) {
  int numSpecies = speciesList.size();
  for (int iSpecies=0; iSpecies<numSpecies; iSpecies++) {
    vector<Potential*> &iBondedPotentials = bondedPotentials[iSpecies];
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
          handleOneBondPair(doForces, uTot, iAtom, jAtom, p);
        }
        firstAtom += sna;
      }
    }

    vector<PotentialAngle*> &iBondAnglePotentials = bondAnglePotentials[iSpecies];
    vector<vector<int*> > iBondAngleTriplets = bondAngleTriplets[iSpecies];
    for (int j=0; j<(int)iBondAnglePotentials.size(); j++) {
      PotentialAngle* p = iBondAnglePotentials[j];
      vector<int*> jBondAngleTriplets = iBondAngleTriplets[j];
      int firstAtom = box.getFirstAtom(iSpecies, 0);
      for (int kMolecule=0; kMolecule<box.getNumMolecules(iSpecies); kMolecule++) {
        for (int l=0; l<(int)jBondAngleTriplets.size(); l++) {
          int* lBondAngleTriplets = jBondAngleTriplets[l];
          int iAtom = firstAtom + lBondAngleTriplets[0];
          int jAtom = firstAtom + lBondAngleTriplets[1];
          int kAtom = firstAtom + lBondAngleTriplets[2];
          handleOneBondAngleTriplet(doForces, uTot, iAtom, jAtom, kAtom, p);
        }
        firstAtom += sna;
      }
    }
  }
}

void PotentialMaster::computeOneMoleculeBonds(const int iSpecies, const int iMolecule, double &u1) {
  vector<Potential*> &iBondedPotentials = bondedPotentials[iSpecies];
  if (iBondedPotentials.size() == 0) return;
  vector<vector<int*> > iBondedPairs = bondedPairs[iSpecies];
  int firstAtom = box.getFirstAtom(iSpecies, iMolecule);
  for (int j=0; j<(int)iBondedPotentials.size(); j++) {
    Potential* p = iBondedPotentials[j];
    vector<int*> jBondedPairs = iBondedPairs[j];
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
  }
}

double PotentialMaster::oldEnergy(int iAtom) {
  double u = 2*uAtom[iAtom];
  if (doSingleTruncationCorrection) {
    u += computeOneTruncationCorrection(iAtom);
  }
  if (embeddingPotentials) {
    u += oldEmbeddingEnergy(iAtom);
  }
  if (doEwald) {
    u += ewald->computeFourierAtom(iAtom, true);
  }
  return u;
}

double PotentialMaster::oldEmbeddingEnergy(int iAtom) {
  // just compute all the embedding energies
  int numAtoms = box.getNumAtoms();
  rhoAtomsChanged.clear();
  rhoAtomsChanged.push_back(iAtom);
  int iType = box.getAtomType(iAtom);
  double u = embedF[iType]->f(rhoSum[iAtom]);
  double *ri = box.getAtomPosition(iAtom);
  double zero[3];
  zero[0] = zero[1] = zero[2] = 0;
  for (int jAtom=0; jAtom<numAtoms; jAtom++) {
    if (jAtom==iAtom) continue;
    double *rj = box.getAtomPosition(jAtom);
    double dr[3];
    for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
    box.nearestImage(dr);
    handleOldEmbedding(zero, dr, zero, jAtom, u, box.getAtomType(jAtom));
  }
  return u;
}

double PotentialMaster::oldMoleculeEnergy(int iMolecule) {
  // only works for rigid molecules, handles monatomic EAM
  int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
  double u = 0;
  for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
    u += 2*uAtom[iAtom];
    if (doSingleTruncationCorrection) {
      u += computeOneTruncationCorrection(iAtom);
    }
    if (embeddingPotentials) {
      u += oldEmbeddingEnergy(iAtom);
    }
    // we'll double count any interactions between this atom and image of
    // another atom in the same molecule (we count it now, we'll count it
    // again for that other atom).  this method computes only "up", so
    // subtracting will give us the right result
    u -= oldIntraMoleculeEnergyLS(iAtom, iLastAtom);
  }
  if (doEwald) {
    u += ewald->oneMoleculeFourierEnergy(iMolecule, true);
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
  if (embeddingPotentials) {
    int numAtomsChanged = rhoAtomsChanged.size();
    for (int i=0; i<numAtomsChanged; i++) {
      int iAtom = rhoAtomsChanged[i];
      drhoSum[iAtom] = 0;
    }
    rhoAtomsChanged.clear();
  }
  if (doEwald) {
    ewald->resetAtomDU();
  }
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
  uAtomsChanged.clear();
  duAtomSingle = duAtomMulti = false;
  if (embeddingPotentials) {
    numAtomsChanged = rhoAtomsChanged.size();
    for (int i=0; i<numAtomsChanged; i++) {
      int iAtom = rhoAtomsChanged[i];
      // by the time we get to processAtom(+1), we have
      // the difference.  just use that
      if (coeff==1) {
        rhoSum[iAtom] += drhoSum[iAtom];
      }
      // we when get called with -1, we'll ignore drhoSum
      drhoSum[iAtom] = 0;
    }
    rhoAtomsChanged.clear();
  }
  if (doEwald) {
    ewald->processAtomU(coeff);
  }
}

void PotentialMaster::computeOne(const int iAtom, double &u1) {
  duAtomSingle = true;
  u1 = 0;
  uAtomsChanged.resize(1);
  duAtom.resize(1);
  uAtomsChanged[0] = iAtom;
  duAtom[0] = 0;
  if (embeddingPotentials && rhoAtomsChanged.size()==0) {
    // if we called oldEnergy for this atom, then it will already be in the list
    rhoAtomsChanged.push_back(iAtom);
  }
  int iMolecule = 0, iFirstAtom = 0, iSpecies = 0;
  if (!pureAtoms && !rigidMolecules) {
    box.getMoleculeInfoAtom(iAtom, iMolecule, iSpecies, iFirstAtom);
  }
  const double *ri = box.getAtomPosition(iAtom);
  computeOneInternal(iAtom, ri, u1, iSpecies, iMolecule, iFirstAtom, true);
  if (doEwald) {
    u1 += ewald->computeFourierAtom(iAtom, false);
  }
  if (doSingleTruncationCorrection) {
    u1 += computeOneTruncationCorrection(iAtom);
  }
}

void PotentialMaster::computeOneInternal(const int iAtom, const double *ri, double &u1, const int iSpecies, const int iMolecule, const int iFirstAtom, const bool onlyAtom) {
  vector<int> *iBondedAtoms = nullptr;
  if (!pureAtoms && !rigidMolecules) {
    iBondedAtoms = &bondedAtoms[iSpecies][iAtom-iFirstAtom];
  }
  int iType = box.getAtomType(iAtom);
  double* iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];
  int numAtoms = box.getNumAtoms();
  double iRhoCutoff = embeddingPotentials ? rhoCutoffs[iType] : 0;
  double zero[3];
  zero[0] = zero[1] = zero[2] = 0;
  Potential* iRhoPotential = embeddingPotentials ? rhoPotentials[iType] : nullptr;
  for (int jAtom=0; jAtom<numAtoms; jAtom++) {
    if (jAtom==iAtom) continue;
    if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
    int jType = box.getAtomType(jAtom);
    Potential* pij = iPotentials[jType];
    if (!pij) continue;
    double *rj = box.getAtomPosition(jAtom);
    double dr[3];
    for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
    box.nearestImage(dr);
    handleComputeOne(pij, zero, dr, zero, iAtom, jAtom, u1, iCutoffs[jType], iRhoCutoff, iRhoPotential, iType, jType, false);
  }
  if (embeddingPotentials) {
    // we just computed new rhoSum[iAtom].  now subtract the old one
    drhoSum[iAtom] -= rhoSum[iAtom];
    u1 += embedF[iType]->f(rhoSum[iAtom] + drhoSum[iAtom]);
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
    computeOneInternal(iAtom, ri, u1, iSpecies, iMolecule, firstAtom, false);
    if (doSingleTruncationCorrection) {
      u1 += computeOneTruncationCorrection(iAtom);
    }
  }
  if (!pureAtoms && !rigidMolecules) {
    computeOneMoleculeBonds(iSpecies, iMolecule, u1);
  }
  if (doEwald) {
    u1 += ewald->oneMoleculeFourierEnergy(iMolecule, false);
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
  if (embeddingPotentials) {
    for (int iAtom=0; iAtom<numAtoms; iAtom++) {
      int iType = box.getAtomType(iAtom);
      uTot += embedF[iType]->f(rhoSum[iAtom]);
    }
  }
  if (doEwald) {
    uTot += ewald->uTotalFromAtoms();
  }
  double virialTot;
  computeAllTruncationCorrection(uTot, virialTot);
  return uTot;
}

void PotentialMaster::updateAtom(int iAtom) {
  box.nearestImage(box.getAtomPosition(iAtom));
}
