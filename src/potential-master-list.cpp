#include <stdio.h>
#include "potential-master.h"
#include "alloc2d.h"

PotentialMasterList::PotentialMasterList(const SpeciesList& sl, Box& box, bool doEmbed, int cellRange, double nRange) : PotentialMasterCell(sl, box, doEmbed, cellRange), nbrRange(nRange), nbrs(nullptr), onlyUpNbrs(true), numAtomNbrsUp(nullptr), numAtomNbrsDn(nullptr), nbrsNumAtoms(0), maxNab(0), nbrBoxOffsets(nullptr), forceReallocNbrs(false), oldAtomPositions(nullptr), safetyFac(0.1) {
  maxR2 = (double*)malloc(numAtomTypes*sizeof(double));
  maxR2Unsafe = (double*)malloc(numAtomTypes*sizeof(double));
}

PotentialMasterList::~PotentialMasterList() {
  free(maxR2);
  free(maxR2Unsafe);
  free(numAtomNbrsUp);
  free(numAtomNbrsDn);
  free2D((void**)nbrs);
  free2D((void**)nbrBoxOffsets);
  free2D((void**)oldAtomPositions);
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
}

void PotentialMasterList::setDoDownNbrs(bool doDown) {
  onlyUpNbrs = !doDown;
}

int PotentialMasterList::checkNbrPair(int iAtom, int jAtom, const bool skipIntra, double *ri, double *rj, double rc2, double minR2, double *jbo) {
  double r2 = 0;
  for (int k=0; k<3; k++) {
    double dr = rj[k]+jbo[k]-ri[k];
    r2 += dr*dr;
  }
  if (r2 > rc2 || (skipIntra && r2 < minR2)) return 0;
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
  double minR2 = 0.5*bs[0];
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
        numAtomNbrsDn[jAtom]++;
      }
    }
  }
}

void PotentialMasterList::computeAll(vector<PotentialCallback*> &callbacks) {
  pairCallbacks.resize(0);
  bool doForces = false;
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if (!embeddingPotentials && (*it)->callPair) pairCallbacks.push_back(*it);
    if ((*it)->takesForces) doForces = true;
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
  double uTot=0, virialTot=0;
  for (int i=0; i<numAtoms; i++) {
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
      handleComputeAll(iAtom, jAtom, ri, rj, jbo, pij, uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, rc2, iRhoPotential, iRhoCutoff, iType, jType, doForces);
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
    computeAllFourier(doForces, uTot);
  }
  computeAllTruncationCorrection(uTot, virialTot);
  if (!pureAtoms && !rigidMolecules) {
    computeAllBonds(doForces, uTot, virialTot);
  }
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}
