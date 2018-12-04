/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <iostream>
#include "ewald.h"
#include "alloc2d.h"
#include "potential-master.h"
#ifdef VALGRIND_CHECKS
#include "valgrind/memcheck.h"
#endif

PotentialMasterCell::PotentialMasterCell(const SpeciesList& sl, Box& box, bool doEmbed, int cRange) : PotentialMaster(sl, box, doEmbed), cellManager(sl, box, cRange), cellRange(cRange), cellNextAtom(cellManager.cellNextAtom), atomCell(cellManager.atomCell), cellLastAtom(cellManager.cellLastAtom), cellOffsets(cellManager.cellOffsets), wrapMap(cellManager.wrapMap), boxOffsets(cellManager.boxOffsets), lsNeeded(false) {
}

PotentialMasterCell::~PotentialMasterCell() {
}

void PotentialMasterCell::setPairPotential(int iType, int jType, Potential *pij) {
  PotentialMaster::setPairPotential(iType, jType, pij);
  cellManager.setRange(getRange());
}

double PotentialMasterCell::getRange() {
  double range = 0;
  for (int iType=0; iType<numAtomTypes; iType++) {
    for (int jType=iType; jType<numAtomTypes; jType++) {
      if (pairPotentials[iType][jType] == nullptr) continue;
      double rc = pairPotentials[iType][jType]->getCutoff();
      if (rc > range) range = rc;
    }
  }
  if (embeddingPotentials) {
    for (int iType=0; iType<numAtomTypes; iType++) {
      if (rhoPotentials[iType] == nullptr) continue;
      double rc = rhoPotentials[iType]->getCutoff();
      if (rc > range) range = rc;
    }
  }
  return range;
}

void PotentialMasterCell::init() {
  PotentialMaster::init();
  cellManager.init();
#ifdef DEBUG
  uAtom[0] = 0;
#endif
}

void PotentialMasterCell::updateVolume() {
  cellManager.init();
#ifdef DEBUG
  uAtom[0] = 0;
#endif
}

int* PotentialMasterCell::getNumCells() {
  return cellManager.getNumCells();
}

void PotentialMasterCell::updateAtom(int iAtom) {
  box.nearestImage(box.getAtomPosition(iAtom));
  cellManager.updateAtom(iAtom);
}

double PotentialMasterCell::oldEmbeddingEnergy(int iAtom) {
  rhoAtomsChanged.clear();
  rhoAtomsChanged.push_back(iAtom);
  int iType = box.getAtomType(iAtom);
  double u = embedF[iType]->f(rhoSum[iAtom]);

  const int iCell = atomCell[iAtom];

  const double *jbo = boxOffsets[iCell];
  const double *ri = box.getAtomPosition(iAtom);
  for (int jAtom = cellLastAtom[iCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
    if (jAtom!=iAtom) handleOldEmbedding(ri, box.getAtomPosition(jAtom), jbo, jAtom, u, box.getAtomType(jAtom));
  }

  for (vector<int>::const_iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
    int jCell = iCell + *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      handleOldEmbedding(ri, box.getAtomPosition(jAtom), jbo, jAtom, u, box.getAtomType(jAtom));
    }
    jCell = iCell - *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      handleOldEmbedding(ri, box.getAtomPosition(jAtom), jbo, jAtom, u, box.getAtomType(jAtom));
    }
  }
  return u;
}

void PotentialMasterCell::computeAll(vector<PotentialCallback*> &callbacks) {
  const double *bs = box.getBoxSize();
  minR2 = 0.5*bs[0];
  for (int k=1; k<3; k++) minR2 = bs[k]<minR2 ? 0.5*bs[k] : minR2;
  minR2 *= minR2;
  pairCallbacks.resize(0);
  bool doForces = false, doPhi = false, doDFDV = false;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if (!embeddingPotentials && (*it)->callPair) pairCallbacks.push_back(*it);
    if ((*it)->takesForces) doForces = true;
    if ((*it)->takesPhi) doPhi = true;
    if ((*it)->takesDFDV) doDFDV = true;
  }
  const int numAtoms = box.getNumAtoms();
  if (doForces && numAtoms > numForceAtoms) {
    force = (double**)malloc2D(numAtoms, 3, sizeof(double));
    if (embeddingPotentials) {
      idf = (double*)realloc(idf, numAtoms*sizeof(double));
    }
    numForceAtoms = numAtoms;
  }
  if (embeddingPotentials && numAtoms > numRhoSumAtoms) {
    drhoSum.resize(numAtoms);
    numRhoSumAtoms = numAtoms;
    rhoSum = (double*)realloc(rhoSum, numAtoms*sizeof(double));
  }
  double uTot=0, virialTot=0;
#ifdef DEBUG
  vector<double> uCheck;
  vector<double> rhoCheck;
  uCheck.resize(box.getNumAtoms());
  if (embeddingPotentials) rhoCheck.resize(box.getNumAtoms());
#endif
  for (int i=0; i<numAtoms; i++) {
#ifdef DEBUG
    uCheck[i] = uAtom[i];
    if (embeddingPotentials) rhoCheck[i] = rhoSum[i];
#endif
    uAtom[i] = 0;
    if (embeddingPotentials) {
      rhoSum[i] = 0;
      drhoSum[i] = 0;
    }
    if (doForces) for (int k=0; k<3; k++) force[i][k] = 0;
  }
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int iMolecule = iAtom;
    vector<int> *iBondedAtoms = nullptr;
    int iSpecies;
    if (!pureAtoms) {
      int iFirstAtom;
      box.getMoleculeInfoAtom(iAtom, iMolecule, iSpecies, iFirstAtom);
      if (!rigidMolecules) {
        iBondedAtoms = &bondedAtoms[iSpecies][iAtom-iFirstAtom];
      }
    }
    const int iType = box.getAtomType(iAtom);
    const double *iCutoffs = pairCutoffs[iType];
    Potential** iPotentials = pairPotentials[iType];
    const double *ri = box.getAtomPosition(iAtom);
    double *fi = doForces ? force[iAtom] : nullptr;
    Potential* iRhoPotential = embeddingPotentials ? rhoPotentials[iType] : nullptr;
    double iRhoCutoff = embeddingPotentials ? rhoCutoffs[iType] : 0;
    int jAtom=iAtom;
    const double *jbo = boxOffsets[atomCell[iAtom]];
#ifdef VALGRIND_CHECKS
    if (VALGRIND_CHECK_MEM_IS_ADDRESSABLE(jbo, 24)) {
      printf("oops 0 not addressable for %d %d\n", iAtom, atomCell[iAtom]);
    }
#endif
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      bool skipIntra = checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms);
      const int jType = box.getAtomType(jAtom);
      Potential* pij = iPotentials[jType];
      if (!pij) continue;
      const double *rj = box.getAtomPosition(jAtom);
      handleComputeAll(iAtom, jAtom, ri, rj, jbo, pij, uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, iCutoffs[jType], iRhoPotential, iRhoCutoff, iType, jType, doForces, skipIntra);
    }
    const int iCell = atomCell[iAtom];
    for (vector<int>::const_iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      jbo = boxOffsets[jCell];
#ifdef VALGRIND_CHECKS
      if (VALGRIND_CHECK_MEM_IS_ADDRESSABLE(jbo, 24)) {
        printf("oops up not addressable for iCell %d offset %d jCell %d\n", iCell, *it, jCell);
        abort();
      }
#endif
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        bool skipIntra = checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms);
        const int jType = box.getAtomType(jAtom);
        Potential* pij = iPotentials[jType];
        if (!pij) continue;
        const double *rj = box.getAtomPosition(jAtom);
        handleComputeAll(iAtom, jAtom, ri, rj, jbo, pij, uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, iCutoffs[jType], iRhoPotential, iRhoCutoff, iType, jType, doForces, skipIntra);
      }
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
      if (doForces) idf[iAtom] = df;
    }
    if (doForces) {
      for (int iAtom=0; iAtom<numAtoms; iAtom++) {
        int iType = box.getAtomType(iAtom);
        double *ri = box.getAtomPosition(iAtom);
        double iRhoCutoff = rhoCutoffs[iType];
        Potential* iRhoPotential = rhoPotentials[iType];
        double df = idf[iAtom];
        int jAtom = iAtom;
        const double *jbo = boxOffsets[atomCell[iAtom]];
        while ((jAtom = cellNextAtom[jAtom]) > -1) {
          int jType = box.getAtomType(jAtom);
          double *rj = box.getAtomPosition(jAtom);
          handleComputeAllEmbed(iAtom, jAtom, iType, jType, ri, rj, jbo, df, virialTot, iRhoPotential, iRhoCutoff, rdrhoIdx);
        }
        const int iCell = atomCell[iAtom];
        for (vector<int>::const_iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
          int jCell = iCell + *it;
          jbo = boxOffsets[jCell];
          jCell = wrapMap[jCell];
          for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
            const int jType = box.getAtomType(jAtom);
            const double *rj = box.getAtomPosition(jAtom);
            handleComputeAllEmbed(iAtom, jAtom, iType, jType, ri, rj, jbo, df, virialTot, iRhoPotential, iRhoCutoff, rdrhoIdx);
          }
        }
      }
      rdrho.clear();
    }
  }
  if (!pureAtoms && !rigidMolecules) {
    computeAllBonds(doForces, uTot);
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
        fprintf(stderr, "PMC uAtomCheck oops %d %f %f %f\n", i, uCheck[i], uAtom[i], uCheck[i]-uAtom[i]);
        oops=true;
      }
    }
    if (embeddingPotentials) {
      for (int i=0; i<numAtoms; i++) {
        if (fabs(rhoCheck[i] - rhoSum[i]) > 1e-7) {
          fprintf(stderr, "PMC rhoSumCheck oops %d %e %e %e\n", i, rhoCheck[i], rhoSum[i], rhoCheck[i]-rhoSum[i]);
          abort();
        }
      }
    }
    if (oops) abort();
  }
#endif
  computeAllTruncationCorrection(uTot, virialTot);
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}

void PotentialMasterCell::computeOneInternal(const int iAtom, const double *ri, double &u1, const int iSpecies, const int iMolecule, const int iFirstAtom, const bool onlyAtom) {
  const int iType = box.getAtomType(iAtom);
  const double *iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];

  const int iCell = atomCell[iAtom];

  vector<int> *iBondedAtoms = nullptr;
  if (!pureAtoms && !rigidMolecules) {
    iBondedAtoms = &bondedAtoms[iSpecies][iAtom-iFirstAtom];
  }
  double iRhoCutoff = embeddingPotentials ? rhoCutoffs[iType] : 0;
  Potential* iRhoPotential = embeddingPotentials ? rhoPotentials[iType] : nullptr;
  const double *jbo = boxOffsets[iCell];
  for (int jAtom = cellLastAtom[iCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
    if (jAtom!=iAtom) {
      if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
      const int jType = box.getAtomType(jAtom);
      Potential* pij = iPotentials[jType];
      if (!pij) continue;
      const double *rj = box.getAtomPosition(jAtom);
      handleComputeOne(pij, ri, rj, jbo, iAtom, jAtom, u1, iCutoffs[jType], iRhoCutoff, iRhoPotential, iType, jType, false);
    }
  }

  for (vector<int>::const_iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
    int jCell = iCell + *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      bool skipIntra = checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms);
      const int jType = box.getAtomType(jAtom);
      Potential* pij = iPotentials[jType];
      if (!pij) continue;
      const double *rj = box.getAtomPosition(jAtom);
      handleComputeOne(pij, ri, rj, jbo, iAtom, jAtom, u1, iCutoffs[jType], iRhoCutoff, iRhoPotential, iType, jType, skipIntra);
    }
    // now down
    jCell = iCell - *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      bool skipIntra = checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms);
      if (skipIntra && !onlyAtom) continue;
      const int jType = box.getAtomType(jAtom);
      Potential* pij = iPotentials[jType];
      if (!pij) continue;
      const double *rj = box.getAtomPosition(jAtom);
      handleComputeOne(pij, ri, rj, jbo, iAtom, jAtom, u1, iCutoffs[jType], iRhoCutoff, iRhoPotential, iType, jType, skipIntra);
    }
  }
  if (embeddingPotentials) {
    // we just computed new rhoSum[iAtom].  now subtract the old one
    drhoSum[iAtom] -= rhoSum[iAtom];
    u1 += embedF[iType]->f(rhoSum[iAtom] + drhoSum[iAtom]);
  }
}

double PotentialMasterCell::oldIntraMoleculeEnergyLS(int iAtom, int iLastAtom) {
  // pretend we're doing computeOne... handleComputeOne expects this stuff to exist
  duAtomSingle = true;
  uAtomsChanged.resize(1);
  duAtom.resize(1);
  uAtomsChanged[0] = iAtom;
  duAtom[0] = 0;

  const int iType = box.getAtomType(iAtom);
  const double *iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];
  double* ri = box.getAtomPosition(iAtom);

  double** rawBoxOffsets = cellManager.rawBoxOffsets;
  int numRawBoxOffsets = cellManager.numRawBoxOffsets;
  double u = 0;
  for (int jAtom=iAtom; jAtom<=iLastAtom; jAtom++) {
    const int jType = box.getAtomType(jAtom);
    Potential* pij = iPotentials[jType];
    if (!pij) continue;
    double rc2 = iCutoffs[jType];
    if (rc2 < minR2) continue;
    const double *rj = box.getAtomPosition(jAtom);
    for (int ijbo=0; ijbo<numRawBoxOffsets; ijbo++) {
      double* jbo = rawBoxOffsets[ijbo];
      handleComputeOne(pij, ri, rj, jbo, iAtom, jAtom, u, rc2, 0, nullptr, iType, jType, true);
    }
  }
  duAtomSingle = false;
  uAtomsChanged.clear();
  duAtom.clear();
  return u;
}

void PotentialMasterCell::removeMolecule(int iSpecies, int iMolecule) {
  // our base class fixs up uAtom by filling in the hole left by iMolecule with 
  // the data from the last molecule of iSpecies, and then shifts all other
  // species down
  PotentialMaster::removeMolecule(iSpecies, iMolecule);
  // we need to do the same with cell data
  cellManager.removeMolecule(iSpecies, iMolecule);
}

void PotentialMasterCell::newMolecule(int iSpecies) {
  // our base class makes room for our atoms in uAtoms and then shifts
  // later species to make room.
  PotentialMaster::newMolecule(iSpecies);
  // we need to repeat that for cell stuff
  cellManager.newMolecule(iSpecies);
}
