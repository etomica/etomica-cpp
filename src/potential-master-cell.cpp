#include <stdio.h>
#include <iostream>
#include "alloc2d.h"
#include "potential-master.h"

PotentialMasterCell::PotentialMasterCell(const SpeciesList& sl, Box& box, int cRange) : PotentialMaster(sl, box), cellManager(sl, box, cRange), cellRange(cRange), cellNextAtom(cellManager.cellNextAtom), atomCell(cellManager.atomCell), cellLastAtom(cellManager.cellLastAtom), cellOffsets(cellManager.cellOffsets), wrapMap(cellManager.wrapMap), boxOffsets(cellManager.boxOffsets), lsNeeded(false) {
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
      double rc = pairPotentials[iType][jType]->getCutoff();
      if (rc > range) range = rc;
    }
  }
  return range;
}

void PotentialMasterCell::init() {
  cellManager.init();
}

int* PotentialMasterCell::getNumCells() {
  return cellManager.getNumCells();
}

void PotentialMasterCell::updateAtom(int iAtom) {
  cellManager.updateAtom(iAtom);
}

void PotentialMasterCell::computeAll(vector<PotentialCallback*> &callbacks) {
  pairCallbacks.resize(0);
  bool doForces = false;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callPair) pairCallbacks.push_back(*it);
    if ((*it)->takesForces) doForces = true;
  }
  if (doForces && !force) {
    force = (double**)malloc2D(box.getNumAtoms(), 3, sizeof(double));
  }
  const int numAtoms = box.getNumAtoms();
  double uTot=0, virialTot=0;
#ifdef DEBUG
  vector<double> uCheck;
  uCheck.resize(box.getNumAtoms());
#endif
  for (int i=0; i<numAtoms; i++) {
#ifdef DEBUG
    uCheck[i] = uAtom[i];
#endif
    uAtom[i] = 0;
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
    int jAtom=iAtom;
    const double *jbo = boxOffsets[atomCell[iAtom]];
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
      const double *rj = box.getAtomPosition(jAtom);
      const int jType = box.getAtomType(jAtom);
      handleComputeAll(iAtom, jAtom, ri, rj, jbo, iPotentials[jType], uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, iCutoffs[jType], doForces);
    }
    const int iCell = atomCell[iAtom];
    for (vector<int>::const_iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      jbo = boxOffsets[jCell];
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
        const int jType = box.getAtomType(jAtom);
        const double *rj = box.getAtomPosition(jAtom);
        handleComputeAll(iAtom, jAtom, ri, rj, jbo, iPotentials[jType], uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, iCutoffs[jType], doForces);
      }
    }
  }
  if (!pureAtoms && !rigidMolecules) {
    computeAllBonds(doForces, uTot, virialTot);
  }
#ifdef DEBUG
  if (uCheck[0]!=0) {
    bool oops = false;
    for (int i=0; i<numAtoms; i++) {
      if (fabs(uCheck[i]-uAtom[i]) > 1e-7) {
        fprintf(stderr, "PMC uAtomCheck oops %d %f %f %f\n", i, uCheck[i], uAtom[i], uCheck[i]-uAtom[i]);
        oops=true;
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

void PotentialMasterCell::computeOneInternal(const int iAtom, const double *ri, double &u1, const int iSpecies, const int iMolecule, const int iFirstAtom) {
  const int iType = box.getAtomType(iAtom);
  const double *iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];

  const int iCell = atomCell[iAtom];

  vector<int> *iBondedAtoms = nullptr;
  if (!pureAtoms && !rigidMolecules) {
    iBondedAtoms = &bondedAtoms[iSpecies][iAtom-iFirstAtom];
  }
  const double *jbo = boxOffsets[iCell];
  for (int jAtom = cellLastAtom[iCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
    if (jAtom!=iAtom) {
      if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
      const int jType = box.getAtomType(jAtom);
      const double *rj = box.getAtomPosition(jAtom);
      handleComputeOne(iPotentials[jType], ri, rj, jbo, iAtom, jAtom, u1, iCutoffs[jType]);
    }
  }

  for (vector<int>::const_iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
    int jCell = iCell + *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
      const double *rj = box.getAtomPosition(jAtom);
      const int jType = box.getAtomType(jAtom);
      handleComputeOne(iPotentials[jType], ri, rj, jbo, iAtom, jAtom, u1, iCutoffs[jType]);
    }
    jCell = iCell - *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      if (checkSkip(jAtom, iSpecies, iMolecule, iBondedAtoms)) continue;
      const double *rj = box.getAtomPosition(jAtom);
      const int jType = box.getAtomType(jAtom);
      handleComputeOne(iPotentials[jType], ri, rj, jbo, iAtom, jAtom, u1, iCutoffs[jType]);
    }
  }
}

void PotentialMasterCell::removeAtom(int iAtom) {
  cellManager.removeAtom(iAtom);
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
