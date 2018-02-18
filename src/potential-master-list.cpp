#include <stdio.h>
#include "potential-master.h"
#include "alloc2d.h"

PotentialMasterList::PotentialMasterList(Potential& p2, Box& box, double pRange, int cellRange, double nbrRange) : PotentialMasterCell(p2, box, nbrRange, cellRange), nbrs(nullptr), onlyUpNbrs(true), numAtomNbrsUp(nullptr), numAtomNbrsDn(nullptr), nbrsNumAtoms(0), maxNab(0), nbrBoxOffsets(nullptr), forceReallocNbrs(false), potentialRange(pRange), oldAtomPositions(nullptr), safetyFac(0.1) { }

PotentialMasterList::~PotentialMasterList() {
  if (numAtomNbrsUp) free(numAtomNbrsUp);
  if (numAtomNbrsDn) free(numAtomNbrsDn);
  if (nbrs) free2D((void**)nbrs);
  if (nbrBoxOffsets) free2D((void**)nbrBoxOffsets);
  if (oldAtomPositions) free2D((void**)oldAtomPositions);
}

void PotentialMasterList::setDoDownNbrs(bool doDown) {
  onlyUpNbrs = !doDown;
}

int PotentialMasterList::checkNbrPair(int iAtom, int jAtom, double *ri, double *rj, double rc2, double *jbo) {
  double r2 = 0;
  double dr[3];
  for (int k=0; k<3; k++) {
    dr[k] = rj[k]-ri[k];
    r2 += dr[k]*dr[k];
  }
  if (r2 > rc2) return 0;
  if (numAtomNbrsUp[iAtom] >= maxNab) return 1;
  nbrs[iAtom][numAtomNbrsUp[iAtom]] = jAtom;
  nbrBoxOffsets[iAtom][numAtomNbrsUp[iAtom]] = jbo;
  numAtomNbrsUp[iAtom]++;
  return 0;
}

void PotentialMasterList::checkUpdateNbrs() {
  int boxNumAtoms = box.getNumAtoms();
  double maxDrUnsafe = (range-potential.getCutoff())*0.5;
#ifdef DEBUG
  double maxR2Unsafe = maxDrUnsafe*maxDrUnsafe;
#endif
  double maxDr = maxDrUnsafe*(1-safetyFac);
  double maxR2 = maxDr*maxDr;
  bool needsUpdate = false;
  for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    double r2 = 0;
    for (int j=0; j<3; j++) {
      double dr = ri[j]-oldAtomPositions[iAtom][j];
      r2 += dr*dr;
    }
    if (r2 > maxR2) {
#ifdef DEBUG
      if (safetyFac>0 && r2 > maxR2Unsafe) {
        fprintf(stderr, "atom %d drifted into unsafe zone before nbr update\n");
      }
      needsUpdate = true;
#else
      reset();
      return;
#endif
    }
  }
  if (needsUpdate) reset();
}

void PotentialMasterList::reset() {
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

  assignCells();

// if an atom has more neighbors than we've allocated space for, then
// maxNab will be increased and execution will return to this point
resetStart:
  if ((!numAtomNbrsUp || moreAtoms) && boxNumAtoms>0)  {
    numAtomNbrsUp = (int*)realloc(numAtomNbrsUp, boxNumAtoms*sizeof(int));
    if (!onlyUpNbrs) numAtomNbrsDn = (int*)realloc(numAtomNbrsDn, boxNumAtoms*sizeof(int));
  }
  // forceReallocNbrs can be used to force reallocation when max # of nbrs is too small
  if ((!nbrs || moreAtoms || forceReallocNbrs) && boxNumAtoms>0)  {
    if (!maxNab) maxNab = 5;
    if (boxNumAtoms > nbrsNumAtoms) nbrsNumAtoms = boxNumAtoms;
    nbrs = (int**)realloc2D((void**)nbrs, nbrsNumAtoms, maxNab, sizeof(int));
    nbrBoxOffsets = (double***)realloc2D((void**)nbrBoxOffsets, nbrsNumAtoms, maxNab, sizeof(double*));
    forceReallocNbrs = false;
  }

  for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) numAtomNbrsUp[iAtom] = 0;

  double rc2 = range*range;
  double rjp[3];
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int tooMuch = 0;
    double *ri = box.getAtomPosition(iAtom);
    int jAtom=iAtom;
    int iCell = atomCell[iAtom];
    double *jbo = boxOffsets[iCell]; // always 0
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      double *rj = box.getAtomPosition(jAtom);
      tooMuch += checkNbrPair(iAtom, jAtom, ri, rj, rc2, jbo);
    }
    for (vector<int>::iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      jbo = boxOffsets[jCell];
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        double *rj = box.getAtomPosition(jAtom);
        for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
        tooMuch += checkNbrPair(iAtom, jAtom, ri, rjp, rc2, jbo);
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
        numAtomNbrsDn[jAtom]++;
      }
    }
  }
}

void PotentialMasterList::computeAll(vector<PotentialCallback*> &callbacks) {
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
  double uTot=0, virialTot=0;
  double rc = potential.getCutoff();
  double rc2 = rc*rc;
  double rjp[3];
  for (int i=0; i<numAtoms; i++) {
    uAtom[i] = 0;
    if (doForces) for (int k=0; k<3; k++) force[i][k] = 0;
  }
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    double *fi = doForces ? force[iAtom] : nullptr;
    int iNumNbrs = numAtomNbrsUp[iAtom];
    int* iNbrs = nbrs[iAtom];
    double** iNbrBoxOffsets = nbrBoxOffsets[iAtom];
    for (int j=0; j<iNumNbrs; j++) {
      int jAtom = iNbrs[j];
      double *rj = box.getAtomPosition(jAtom);
      double *jbo = iNbrBoxOffsets[j];
      for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
      handleComputeAll(iAtom, jAtom, ri, rjp, uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, rc2, doForces);
    }
  }
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}
