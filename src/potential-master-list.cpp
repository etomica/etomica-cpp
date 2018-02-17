#include "potential-master.h"
#include "alloc2d.h"

PotentialMasterList::PotentialMasterList(Potential& p2, Box& box, double potentialRange, int cellRange, double nbrRange) : PotentialMasterCell(p2, box, potentialRange, cellRange), nbrs(nullptr), onlyUpNbrs(true), numAtomNbrsUp(nullptr), numAtomNbrsDn(nullptr), nbrsNumAtoms(0), maxNab(0), forceReallocNbrs(0) {
}

PotentialMasterList::~PotentialMasterList() {
  if (numAtomNbrsUp) free(numAtomNbrsUp);
  if (numAtomNbrsDn) free(numAtomNbrsDn);
  if (nbrs) free2D((void**)nbrs);
}

void PotentialMasterList::checkNbrPair(int iAtom, int jAtom, double *ri, double *rj, double rc2) {
  double r2 = 0;
  double dr[3];
  for (int k=0; k<3; k++) {
    dr[k] = rj[k]-ri[k];
    r2 += dr[k]*dr[k];
  }
  if (r2 > rc2) return;
  nbrs[iAtom][numAtomNbrsUp[iAtom]] = jAtom;
  numAtomNbrsUp[iAtom]++;
}

void PotentialMasterList::reset() {
  int boxNumAtoms = box.getNumAtoms();
  if ((!numAtomNbrsUp || boxNumAtoms > nbrsNumAtoms) && boxNumAtoms>0)  {
    numAtomNbrsUp = (int*)realloc(numAtomNbrsUp, boxNumAtoms*sizeof(int));
    if (!onlyUpNbrs) numAtomNbrsDn = (int*)realloc(numAtomNbrsDn, boxNumAtoms*sizeof(int));
  }
  // forceReallocNbrs can be used to force reallocation when max # of nbrs is too small
  if ((!nbrs || boxNumAtoms > nbrsNumAtoms || forceReallocNbrs) && boxNumAtoms>0)  {
    if (!maxNab) maxNab = 5;
    if (boxNumAtoms > nbrsNumAtoms) nbrsNumAtoms = boxNumAtoms;
    nbrs = (int**)realloc2D((void**)nbrs, nbrsNumAtoms, maxNab, sizeof(int));
    forceReallocNbrs = false;
  }

  assignCells();

  for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) numAtomNbrsUp[iAtom] = 0;
  if (!onlyUpNbrs) for (int iAtom=0; iAtom<boxNumAtoms; iAtom++) numAtomNbrsDn[iAtom] = 0;

  double rc = potential.getCutoff();
  double rc2 = rc*rc;
  double rjp[3];
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    double *ri = box.getAtomPosition(iAtom);
    int jAtom=iAtom;
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      double *rj = box.getAtomPosition(jAtom);
      checkNbrPair(iAtom, jAtom, ri, rj, rc2);
    }
    int iCell = atomCell[iAtom];
    for (vector<int>::iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      double *jbo = boxOffsets[jCell];
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        double *rj = box.getAtomPosition(jAtom);
        for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
        checkNbrPair(iAtom, jAtom, ri, rjp, rc2);
      }
    }
    if (numAtomNbrsUp[iAtom] > maxNab) {
      maxNab = numAtomNbrsUp[iAtom];
      forceReallocNbrs = true;
      reset();
      return;
    }
  }
}
