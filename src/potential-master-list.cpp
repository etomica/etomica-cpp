/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include "potential-master.h"
#include "ewald.h"
#include "alloc2d.h"

PotentialMasterList::PotentialMasterList(const SpeciesList& sl, Box& box, bool doEmbed, int cellRange, double nRange) : PotentialMasterCell(sl, box, doEmbed, cellRange), nbrRange(nRange), nbrs(nullptr), onlyUpNbrs(true), numAtomNbrsUp(nullptr), numAtomNbrsDn(nullptr), nbrsNumAtoms(0), maxNab(0), nbrBoxOffsets(nullptr), oldAtomPositions(nullptr), safetyFac(0.1) {
  maxR2 = (double*)malloc(numAtomTypes*sizeof(double));
  maxR2Unsafe = (double*)malloc(numAtomTypes*sizeof(double));
  moleculeCutoffs = nullptr;
}

PotentialMasterList::~PotentialMasterList() {
  free(maxR2);
  free(maxR2Unsafe);
  free(numAtomNbrsUp);
  free(numAtomNbrsDn);
  free2D((void**)nbrs);
  free2D((void**)nbrBoxOffsets);
  free2D((void**)oldAtomPositions);
  free2D((void**)moleculeCutoffs);
}

double PotentialMasterList::getRange() {
  return nbrRange;
}

void PotentialMasterList::init() {
  PotentialMasterCell::init();
  double maxRhoCut = 0;
  if (embeddingPotentials) {
    for (int i=0; i<numAtomTypes; i++) {
      if (rhoCutoffs[i] > maxRhoCut) maxRhoCut = rhoCutoffs[i];
    }
  }
  for (int i=0; i<numAtomTypes; i++) {
    maxR2Unsafe[i] = maxR2[i] = 1e100;
    for (int j=0; j<numAtomTypes; j++) {
      double rc = sqrt(pairCutoffs[i][j]);
      if (maxRhoCut > rc) rc = maxRhoCut;
      double maxDrUnsafe = (nbrRange-rc)*0.5;
      double x = maxDrUnsafe*maxDrUnsafe;
      if (maxR2Unsafe[i] < x) continue;
      maxR2Unsafe[i] = x;
      double maxDr = maxDrUnsafe*(1-safetyFac);
      maxR2[i] = maxDr*maxDr;
    }
  }
  reset();
}

void PotentialMasterList::setDoDownNbrs(bool doDown) {
  onlyUpNbrs = !doDown;
}

int PotentialMasterList::checkNbrPair(int iAtom, int jAtom, const bool skipIntra, double *ri, double *rj, double rc2, double minR2, double *jbo) {
  if (moleculeCutoffs) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom;
    box.getMoleculeInfoAtom(iAtom, iMoleculeInSpecies, iSpecies, iFirstAtom);
    int iMolecule = box.getGlobalMoleculeIndex(iSpecies, iMoleculeInSpecies);
    int jSpecies, jMoleculeInSpecies, jFirstAtom;
    box.getMoleculeInfoAtom(jAtom, jMoleculeInSpecies, jSpecies, jFirstAtom);
    int jMolecule = box.getGlobalMoleculeIndex(jSpecies, jMoleculeInSpecies);
    if (iMolecule != jMolecule) {

      if (moleculeNotNbrs[iMolecule].count(jMolecule) > 0) return 0;

      if (moleculeNbrs[iMolecule].count(jMolecule) == 0) {
        Species* is = speciesList.get(iSpecies);
        Species* js = speciesList.get(jSpecies);
        double dr[3];
        double* c = is->getMoleculeCOM(box, iFirstAtom, iFirstAtom+is->getNumAtoms()-1);
        for (int k=0; k<3; k++) dr[k] = c[k];
        c = js->getMoleculeCOM(box, jFirstAtom, jFirstAtom+js->getNumAtoms()-1);

        for (int k=0; k<3; k++) dr[k] -= c[k];
        // incompatible with a lattice sum
        // need to handle not only center being on opposite side from atom (and so jbo being wrong)
        // but also need to remember that i and j are neighbors only for this jbo
        box.nearestImage(dr);
        double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        if (r2 > moleculeCutoffs[iSpecies][jSpecies]*moleculeCutoffs[iSpecies][jSpecies]) {
          moleculeNotNbrs[iMolecule].insert(jMolecule);
          return 0;
        }
        moleculeNbrs[iMolecule].insert(jMolecule);
      }
    }
  }
  else {
    double r2 = 0;
    for (int k=0; k<3; k++) {
      double dr = rj[k]+jbo[k]-ri[k];
      r2 += dr*dr;
    }
    if (r2 > rc2 || (skipIntra && r2 < minR2)) return 0;
  }
  if (numAtomNbrsUp[iAtom] >= maxNab) return 1;
  nbrs[iAtom][numAtomNbrsUp[iAtom]] = jAtom;
  nbrBoxOffsets[iAtom][numAtomNbrsUp[iAtom]] = jbo;
  numAtomNbrsUp[iAtom]++;
  return 0;
}

void PotentialMasterList::checkUpdateNbrs() {
  int boxNumAtoms = box.getNumAtoms();
#ifdef DEBUG
  bool needsUpdate = false;
#endif
  for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    double r2 = 0;
    int iType = box.getAtomType(iAtom);
    for (int j=0; j<3; j++) {
      double dr = ri[j]-oldAtomPositions[iAtom][j];
      r2 += dr*dr;
    }
    if (r2 > maxR2[iType]) {
#ifdef DEBUG
      if (safetyFac>0 && r2 > maxR2Unsafe[iType]) {
        fprintf(stderr, "atom %d drifted into unsafe zone before nbr update\n", iAtom);
      }
      needsUpdate = true;
#else
      reset();
      return;
#endif
    }
  }
#ifdef DEBUG
  if (needsUpdate) reset();
#endif
}

void PotentialMasterList::reset() {
  if (moleculeCutoffs) {
    int nm = box.getTotalNumMolecules();
    moleculeNbrs.resize(nm);
    moleculeNotNbrs.resize(nm);
    for (int i=0; i<nm; i++) {
      moleculeNbrs[i].clear();
      moleculeNotNbrs[i].clear();
    }
  }
  int boxNumAtoms = box.getNumAtoms();
  if (boxNumAtoms==0) return;
  bool moreAtoms = boxNumAtoms > nbrsNumAtoms;
  if (moreAtoms) {
    oldAtomPositions = (double**)realloc2D((void**)oldAtomPositions, boxNumAtoms, 3, sizeof(double));
  }
  for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    box.nearestImage(ri);
    for (int j=0; j<3; j++) oldAtomPositions[iAtom][j] = ri[j];
  }

  cellManager.assignCells();

// if an atom has more neighbors than we've allocated space for, then
// maxNab will be increased and execution will return to this point
  bool forceReallocNbrs = false;
resetStart:
  if ((!numAtomNbrsUp || moreAtoms) && boxNumAtoms>0)  {
    numAtomNbrsUp = (int*)realloc(numAtomNbrsUp, boxNumAtoms*sizeof(int));
    if (!onlyUpNbrs) numAtomNbrsDn = (int*)realloc(numAtomNbrsDn, boxNumAtoms*sizeof(int));
  }
  // forceReallocNbrs can be used to force reallocation when max # of nbrs is too small
  if ((!nbrs || moreAtoms || forceReallocNbrs) && boxNumAtoms>0)  {
    if (!maxNab) maxNab = 5;
    if (moreAtoms) nbrsNumAtoms = boxNumAtoms;
    nbrs = (int**)realloc2D((void**)nbrs, nbrsNumAtoms, maxNab, sizeof(int));
    nbrBoxOffsets = (double***)realloc2D((void**)nbrBoxOffsets, nbrsNumAtoms, maxNab, sizeof(double*));
    forceReallocNbrs = false;
  }

  for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) numAtomNbrsUp[iAtom] = 0;

  double rc2 = nbrRange*nbrRange;
  const double *bs = box.getBoxSize();
  minR2 = 0.5*bs[0];
  for (int k=1; k<3; k++) minR2 = bs[k]<minR2 ? 0.5*bs[k] : minR2;
  minR2 *= minR2;
  for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) {
    int iMolecule = iAtom;
    vector<int> *iBondedAtoms = nullptr;
    int iSpecies = 0;
    if (!pureAtoms) {
      int iFirstAtom;
      box.getMoleculeInfoAtom(iAtom, iMolecule, iSpecies, iFirstAtom);
      if (!rigidMolecules) {
        iBondedAtoms = &bondedAtoms[iSpecies][iAtom-iFirstAtom];
      }
    }
    int tooMuch = 0;
    double *ri = box.getAtomPosition(iAtom);
    int jAtom=iAtom;
    int iCell = atomCell[iAtom];
    double *jbo = boxOffsets[iCell]; // always 0
    Potential** iPotentials = pairPotentials[box.getAtomType(iAtom)];
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      if (!iPotentials[box.getAtomType(jAtom)]) continue;
      if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
      double *rj = box.getAtomPosition(jAtom);
      tooMuch += checkNbrPair(iAtom, jAtom, false, ri, rj, rc2, minR2, jbo);
    }
    for (vector<int>::const_iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      jbo = boxOffsets[jCell];
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        if (!iPotentials[box.getAtomType(jAtom)]) {
          continue;
        }
        bool skipIntra = checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms);
        double *rj = box.getAtomPosition(jAtom);
        tooMuch += checkNbrPair(iAtom, jAtom, skipIntra, ri, rj, rc2, minR2, jbo);
      }
    }
    if (tooMuch > 0) {
#ifdef DEBUG
      printf("maxNab %d => %d\n", maxNab, maxNab+tooMuch);
#endif
      maxNab += tooMuch;
      forceReallocNbrs = true;
      goto resetStart;
    }
  }
  if (!onlyUpNbrs) {
    for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) numAtomNbrsDn[iAtom] = 0;
    for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) {
      int iNumNbrs = numAtomNbrsUp[iAtom];
      int* iNbrs = nbrs[iAtom];
      for (int j=0; j<iNumNbrs; j++) {
        int jAtom = iNbrs[j];
        if (numAtomNbrsDn[jAtom]+numAtomNbrsUp[jAtom] >= maxNab) {
#ifdef DEBUG
          printf("maxNab %d => %d\n", maxNab, maxNab*4/3);
#endif
          maxNab = (maxNab*4) / 3;
          forceReallocNbrs = true;
          goto resetStart;
        }

        nbrs[jAtom][maxNab-1-numAtomNbrsDn[jAtom]] = iAtom;
        // these offsets are actually opposite... we don't have a clean
        // way to find the opposite offset.
        nbrBoxOffsets[jAtom][maxNab-1-numAtomNbrsDn[jAtom]] = nbrBoxOffsets[iAtom][j];
        numAtomNbrsDn[jAtom]++;
      }
    }
  }
}

void PotentialMasterList::computeAll(vector<PotentialCallback*> &callbacks) {
  pairCallbacks.resize(0);
  bool doForces = false, doPhi = false, doDFDV = false;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if (!embeddingPotentials && (*it)->callPair) pairCallbacks.push_back(*it);
    if ((*it)->takesForces) doForces = true;
    if ((*it)->takesPhi) doPhi = true;
    if ((*it)->takesDFDV) doDFDV = true;
  }
  int numAtoms = box.getNumAtoms();
  if (doForces && numAtoms > numForceAtoms) {
    force = (double**)realloc2D((void**)force, numAtoms, 3, sizeof(double));
    if (embeddingPotentials) {
      rhoSum = (double*)realloc(rhoSum, numAtoms*sizeof(double));
      idf = (double*)realloc(idf, numAtoms*sizeof(double));
    }
    numForceAtoms = numAtoms;
  }
#ifdef DEBUG
  vector<double> uCheck;
  vector<double> rhoCheck;
  uCheck.resize(box.getNumAtoms());
  if (embeddingPotentials) rhoCheck.resize(box.getNumAtoms());
#endif
  double uTot=0, virialTot=0;
  for (int i=0; i<numAtoms; i++) {
#ifdef DEBUG
    uCheck[i] = uAtom[i];
    if (embeddingPotentials) rhoCheck[i] = rhoSum[i];
#endif
    uAtom[i] = 0;
    if (embeddingPotentials) rhoSum[i] = 0;
    if (doForces) for (int k=0; k<3; k++) force[i][k] = 0;
  }
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    int iType = box.getAtomType(iAtom);
    double *iCutoffs = pairCutoffs[iType];
    Potential** iPotentials = pairPotentials[iType];
    double *fi = doForces ? force[iAtom] : nullptr;
    int iNumNbrs = numAtomNbrsUp[iAtom];
    int* iNbrs = nbrs[iAtom];
    double** iNbrBoxOffsets = nbrBoxOffsets[iAtom];
    Potential* iRhoPotential = embeddingPotentials ? rhoPotentials[iType] : nullptr;
    double iRhoCutoff = embeddingPotentials ? rhoCutoffs[iType] : 0;
    for (int j=0; j<iNumNbrs; j++) {
      int jAtom = iNbrs[j];
      int jType = box.getAtomType(jAtom);
      double rc2 = iCutoffs[jType];
      Potential* pij = iPotentials[jType];
      double *rj = box.getAtomPosition(jAtom);
      double *jbo = iNbrBoxOffsets[j];
      handleComputeAll(iAtom, jAtom, ri, rj, jbo, pij, uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, rc2, iRhoPotential, iRhoCutoff, iType, jType, doForces, false);
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
        idf[iAtom] = df;
      }
    }
    if (doForces) {
      for (int iAtom=0; iAtom<numAtoms; iAtom++) {
        int iType = box.getAtomType(iAtom);
        if (doForces) {
          double *ri = box.getAtomPosition(iAtom);
          double iRhoCutoff = rhoCutoffs[iType];
          Potential* iRhoPotential = rhoPotentials[iType];
          int iNumNbrs = numAtomNbrsUp[iAtom];
          int* iNbrs = nbrs[iAtom];
          double** iNbrBoxOffsets = nbrBoxOffsets[iAtom];
          double df = idf[iAtom];
          for (int j=0; j<iNumNbrs; j++) {
            int jAtom = iNbrs[j];
            int jType = box.getAtomType(jAtom);
            double *rj = box.getAtomPosition(jAtom);
            double *jbo = iNbrBoxOffsets[j];

            handleComputeAllEmbed(iAtom, jAtom, iType, jType, ri, rj, jbo, df, virialTot, iRhoPotential, iRhoCutoff, rdrhoIdx);
          }
        }
      }
      rdrho.clear();
    }
  }
  if (doEwald) {
    ewald->computeAllFourier(doForces, doPhi, doDFDV, uTot, virialTot, force, &pairCallbacks);
  }
  if (doForces && !pureAtoms) {
    virialTot += computeVirialIntramolecular();
  }
#ifdef DEBUG
  if (uCheck[0]!=0) {
    bool oops = false;
    for (int i=0; i<numAtoms; i++) {
      if (fabs(uCheck[i]-uAtom[i]) > 1e-7) {
        // we've recomputed the energy for this atom and our previous estimate was wrong
        fprintf(stderr, "PML uAtomCheck oops %d %f %f %f\n", i, uCheck[i], uAtom[i], uCheck[i]-uAtom[i]);
        oops=true;
      }
    }
    if (embeddingPotentials) {
      for (int i=0; i<numAtoms; i++) {
        if (fabs(rhoCheck[i] - rhoSum[i]) > 1e-7) {
          fprintf(stderr, "PML rhoSumCheck oops %d %e %e %e\n", i, rhoCheck[i], rhoSum[i], rhoCheck[i]-rhoSum[i]);
          abort();
        }
      }
    }
    if (oops) abort();
  }
#endif
  computeAllTruncationCorrection(uTot, virialTot);
  if (!pureAtoms && !rigidMolecules) {
    computeAllBonds(doForces, uTot);
  }
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}

// do nothing.  not going to update cell or wrap inside the box.
// both of those things will happen during reset.  doing it before
// that will just confuse us.
void PotentialMasterList::updateAtom(int iAtom) {
}

void PotentialMasterList::computeOneInternal(const int iAtom, const double *ri, double &u1, const int iSpecies, const int iMolecule, const int iFirstAtom, const bool onlyAtom) {
  if (onlyUpNbrs) {
    fprintf(stderr, "down neighbors must be enabled to use computeOne\n");
    abort();
  }
  const int iType = box.getAtomType(iAtom);
  const double *iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];

  double iRhoCutoff = embeddingPotentials ? rhoCutoffs[iType] : 0;
  Potential* iRhoPotential = embeddingPotentials ? rhoPotentials[iType] : nullptr;
  int iNumNbrs = numAtomNbrsUp[iAtom];
  int* iNbrs = nbrs[iAtom];
  double** iNbrBoxOffsets = nbrBoxOffsets[iAtom];
  for (int j=0; j<iNumNbrs; j++) {
    int jAtom = iNbrs[j];
    int jType = box.getAtomType(jAtom);
    Potential* pij = iPotentials[jType];
    double *rj = box.getAtomPosition(jAtom);
    double *jbo = iNbrBoxOffsets[j];
    handleComputeOne(pij, ri, rj, jbo, iAtom, jAtom, u1, iCutoffs[jType], iRhoCutoff, iRhoPotential, iType, jType, false);
  }
  iNumNbrs = numAtomNbrsDn[iAtom];
        //nbrs[jAtom][maxNab-1-numAtomNbrsDn[jAtom]] = iAtom;
  for (int j=maxNab-1; j>maxNab-1-iNumNbrs; j--) {
    int jAtom = iNbrs[j];
    int jType = box.getAtomType(jAtom);
    Potential* pij = iPotentials[jType];
    double *rj = box.getAtomPosition(jAtom);
    double *jbo = iNbrBoxOffsets[j];
    double ibo[3] = {-jbo[0], -jbo[1], -jbo[2]};
    handleComputeOne(pij, ri, rj, ibo, iAtom, jAtom, u1, iCutoffs[jType], iRhoCutoff, iRhoPotential, iType, jType, false);
  }
  if (embeddingPotentials) {
    // we just computed new rhoSum[iAtom].  now subtract the old one
    drhoSum[iAtom] -= rhoSum[iAtom];
    u1 += embedF[iType]->f(rhoSum[iAtom] + drhoSum[iAtom]);
  }
}

void PotentialMasterList::setMoleculePair(int iSpecies, int jSpecies, double rc) {
  if (!moleculeCutoffs) {
    moleculeCutoffs = (double**)malloc2D(numAtomTypes, numAtomTypes, sizeof(double));
    for (int i=0; i<numAtomTypes; i++) {
      for (int j=0; j<numAtomTypes; j++) {
        moleculeCutoffs[i][j] = 0;
      }
    }
  }
  moleculeCutoffs[iSpecies][jSpecies] = moleculeCutoffs[jSpecies][iSpecies] = rc;
}
