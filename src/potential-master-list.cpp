#include <stdio.h>
#include "potential-master.h"
#include "alloc2d.h"

PotentialMasterList::PotentialMasterList(Potential& p2, Box& box, double pRange, int cellRange, double nbrRange) : PotentialMasterCell(p2, box, nbrRange, cellRange), nbrs(nullptr), onlyUpNbrs(true), numAtomNbrsUp(nullptr), numAtomNbrsDn(nullptr), nbrsNumAtoms(0), maxNab(0), forceReallocNbrs(false), potentialRange(pRange) { }

PotentialMasterList::~PotentialMasterList() {
  if (numAtomNbrsUp) free(numAtomNbrsUp);
  if (numAtomNbrsDn) free(numAtomNbrsDn);
  if (nbrs) free2D((void**)nbrs);
}

void PotentialMasterList::setDoDownNbrs(bool doDown) {
  onlyUpNbrs = !doDown;
}

int PotentialMasterList::checkNbrPair(int iAtom, int jAtom, double *ri, double *rj, double rc2) {
  double r2 = 0;
  double dr[3];
  for (int k=0; k<3; k++) {
    dr[k] = rj[k]-ri[k];
    r2 += dr[k]*dr[k];
  }
  if (r2 > rc2) return 0;
  if (numAtomNbrsUp[iAtom] >= maxNab) return 1;
  nbrs[iAtom][numAtomNbrsUp[iAtom]] = jAtom;
  numAtomNbrsUp[iAtom]++;
  return 0;
}

void PotentialMasterList::reset() {
  int boxNumAtoms = box.getNumAtoms();
// if an atom has more neighbors than we've allocated space for, then
// maxNab will be increased and execution will return to this point
resetStart:
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

  double rc2 = range*range;
  double rjp[3];
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int tooMuch = 0;
    double *ri = box.getAtomPosition(iAtom);
    int jAtom=iAtom;
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      double *rj = box.getAtomPosition(jAtom);
      tooMuch += checkNbrPair(iAtom, jAtom, ri, rj, rc2);
    }
    int iCell = atomCell[iAtom];
    for (vector<int>::iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      double *jbo = boxOffsets[jCell];
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        double *rj = box.getAtomPosition(jAtom);
        for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
        tooMuch += checkNbrPair(iAtom, jAtom, ri, rjp, rc2);
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
