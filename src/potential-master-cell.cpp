#include <stdio.h>
#include <iostream>
#include "alloc2d.h"
#include "potential-master.h"

PotentialMasterCell::PotentialMasterCell(const SpeciesList& sl, Box& box, int cRange) : PotentialMaster(sl, box), cellRange(cRange), rawBoxOffsets(nullptr), boxOffsets(nullptr) {
}

PotentialMasterCell::~PotentialMasterCell() {
  if (rawBoxOffsets) free2D((void**)rawBoxOffsets);
  if (boxOffsets) free(boxOffsets);
}

#define i_cell(ax,ay,az) ((az % numCells[2]) \
                        + (ay % numCells[1]) * numCells[2] \
                        + (ax % numCells[0]) * numCells[1] * numCells[2])

#define dx (o[ix])
#define dy (o[iy])
#define dz (oz[iz])

int PotentialMasterCell::wrappedIndex(int i, int nc) {
  int rv = i;
  if (i < cellRange) {
    rv += nc - 2*cellRange;
  }
  else if (nc - i <= cellRange) {
    rv -= nc - 2*cellRange;
  }
  return rv;
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
  double range = getRange();
  double minCellSize = range/cellRange;
  int totalCells = 1;
  double* bs = box.getBoxSize();
  for (int i=0; i<3; i++) {
    // with cell lists, we can accommodate rc>L/2
    // we need the box to be at least the size of a cell.
    // when rc>L/2, we'll end up doing a lattice sum
    if (bs[i] < minCellSize) {
      fprintf(stderr, "box not big enough to accomodate even 1 cell (%f < %f)\n", bs[i], minCellSize);
      exit(1);
    }
    // include cellRange of padding on each side of the box
    numCells[i] = ((int)floor(bs[i]/minCellSize)) + cellRange*2;
    totalCells *= numCells[i];
  }
  cellLastAtom.resize(totalCells);
  wrapMap.resize(totalCells);
  boxOffsets = (double**)realloc(boxOffsets, totalCells*sizeof(double*));
  rawBoxOffsets = (double**)realloc2D((void**)rawBoxOffsets, 3*3*3, 3, sizeof(double));

  cellOffsets.resize(0);
  int dCell = 0;
  int lastCellCount = 0;

  for (int icd=1; icd<=cellRange+5; icd++) {
    // we want all cells whose squared (index) distance is not more than icd
    // but exclude any cells handle in the previous passes
    int icd2 = icd*icd;
    int iz2;
    for (int iz=0; (iz2=iz*iz)<=icd2; iz++) {
      int izm1 = iz==0 ? 0 : (abs(iz)-1);
      int izm1Sq = izm1*izm1;
      int iy2Max = icd2-iz2;
      int iyMax = (int)sqrt(iy2Max+0.001);
      for (int iy=-iyMax; iy<=iyMax; iy++) {
        int iym1 = iy==0 ? 0 : (abs(iy)-1);
        int iym1Sq = iym1*iym1;
        int iy2 = iy*iy;
        int ix2Max = iy2Max-iy2;
        int ixMax = (int)sqrt(ix2Max+0.001);
        int ix2Min = (icd-1)*(icd-1) - iz2 - iy2;
        int ixMin;
        if (ix2Min<0) {
          ix2Min=0;
          ixMin=-1;
        }
        else {
          ixMin = (int)sqrt(ix2Min+0.001);
        }
        for (int ix=-ixMax; ix<=ixMax; ix++) {
          if (ix>=-ixMin && ix<=ixMin) {ix = ixMin+1; if (ix>ixMax) break;}
          if (iz==0) {
            if (iy<1) {if (ix<1) continue;}
            else if (ix<0) continue;
          }
          int ixm1 = ix==0 ? 0 : (abs(ix)-1);
          int ixm1Sq = ixm1*ixm1;
          if (ixm1Sq+iym1Sq+izm1Sq >= cellRange*cellRange) continue;
          int mv = ix+(iy+iz*numCells[1])*numCells[0];
          cellOffsets.push_back(mv);
          dCell++;
        }
      }
    }
    if (dCell==lastCellCount) break;
    lastCellCount = dCell;
  }

  for (int ix=-1; ix<=1; ix++) {
    for (int iy=-1; iy<=1; iy++) {
      for (int iz=-1; iz<=1; iz++) {
        int idx = (ix+1)*9+(iy+1)*3+(iz+1);
        rawBoxOffsets[idx][0] = ix*bs[0];
        rawBoxOffsets[idx][1] = iy*bs[1];
        rawBoxOffsets[idx][2] = iz*bs[2];
      }
    }
  }

  for (int ix=0; ix<numCells[0]; ix++) {
    int x2 = wrappedIndex(ix, numCells[0]);
    int xbo = (ix-x2)/(numCells[0]-2*cellRange);
    for (int iy=0; iy<numCells[1]; iy++) {
      int y2 = wrappedIndex(iy, numCells[1]);
      int ybo = (iy-y2)/(numCells[1]-2*cellRange);
      for (int iz=0; iz<numCells[1]; iz++) {
        int z2 = wrappedIndex(iz, numCells[2]);
        int zbo = (iz-z2)/(numCells[2]-2*cellRange);
        int iMap = i_cell(ix, iy, iz);
        int dCell = i_cell(x2, y2, z2);
        wrapMap[iMap] = dCell;
        boxOffsets[iMap] = rawBoxOffsets[(xbo+1)*9+(ybo+1)*3+(zbo+1)];
      }
    }
  }

  assignCells();
}

int* PotentialMasterCell::getNumCells() {
  return numCells;
}

int PotentialMasterCell::cellForCoord(const double *r) {
  int cellNum = 0;
  const double* bs = box.getBoxSize();
  for (int i=0; i<3; i++) {
    const double x = (r[i] + boxHalf[i])/bs[i];
    cellNum += ((int)(cellRange + x*(numCells[i]-2*cellRange)))*jump[i];
  }
  return cellNum;
}

void PotentialMasterCell::assignCells() {
  double *bs = box.getBoxSize();
  if (bs[0]*bs[1]*bs[2] == 0) {
    fprintf(stderr, "box has 0 volume, can't assign cells\n");
    return;
  }
  if (getRange() == 0 || cellRange == 0) {
    fprintf(stderr, "range and cell range need to be non-zero\n");
    return;
  }

  numAtoms = box.getNumAtoms();
  cellNextAtom.resize(numAtoms);
  atomCell.resize(numAtoms);
  for (int i=0; i<3; i++) boxHalf[i] = 0.5*bs[i];
  jump[0] = numCells[1]*numCells[2];
  jump[1] = numCells[2];
  jump[2] = 1;
  fill(cellLastAtom.begin(), cellLastAtom.end(), -1);
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int cellNum = 0;
    double *r = box.getAtomPosition(iAtom);
    for (int i=0; i<3; i++) {
      double x = (r[i] + boxHalf[i])/bs[i];
      cellNum += ((int)(cellRange + x*(numCells[i]-2*cellRange)))*jump[i];
    }
    atomCell[iAtom] = cellNum;
    cellNextAtom[iAtom] = cellLastAtom[cellNum];
    cellLastAtom[cellNum] = iAtom;
  }
}

void PotentialMasterCell::updateAtom(int iAtom) {
  int cellNum = 0;
  double *bs = box.getBoxSize();
  double *r = box.getAtomPosition(iAtom);
  for (int i=0; i<3; i++) {
    double x = (r[i] + boxHalf[i])/bs[i];
    cellNum += ((int)(cellRange + x*(numCells[i]-2*cellRange)))*jump[i];
  }
  if (iAtom >= numAtoms) {
    int oldNumAtoms = numAtoms;
    numAtoms = box.getNumAtoms();
    //printf("numAtoms => %d\n", numAtoms);
    cellNextAtom.resize(numAtoms);
    atomCell.resize(numAtoms);
    fill(atomCell.begin()+oldNumAtoms, atomCell.end(), -1);
    fill(cellNextAtom.begin()+oldNumAtoms, cellNextAtom.end(), -1);
  }

  int oldCell = atomCell[iAtom];
  // check if existing assignment is right
  if (cellNum == oldCell) return;
  if (oldCell>-1) {
    // delete from old cell
    int j = cellLastAtom[oldCell];
    if (j == iAtom) {
      cellLastAtom[oldCell] = cellNextAtom[iAtom];
    }
    else {
      while (cellNextAtom[j] != iAtom) {
        j = cellNextAtom[j];
      }
      cellNextAtom[j] = cellNextAtom[iAtom];
    }
  }

  atomCell[iAtom] = cellNum;
  cellNextAtom[iAtom] = cellLastAtom[cellNum];
  cellLastAtom[cellNum] = iAtom;
}

void PotentialMasterCell::handleComputeAll(int iAtom, int jAtom, const double *ri, const double *rj, Potential* pij, double &ui, double &uj, double* fi, double* fj, double& uTot, double& virialTot, double rc2, bool doForces) {
  double r2 = 0;
  double dr[3];
  for (int k=0; k<3; k++) {
    dr[k] = rj[k]-ri[k];
    r2 += dr[k]*dr[k];
  }
  if (r2 > rc2) return;
  double u, du, d2u;
  pij->u012(r2, u, du, d2u);
  ui += 0.5*u;
  uj += 0.5*u;
  uTot += u;
  virialTot += du;
  for (vector<PotentialCallback*>::iterator it = pairCallbacks.begin(); it!=pairCallbacks.end(); it++) {
    (*it)->pairCompute(iAtom, jAtom, dr, u, du, d2u);
  }

  // f0 = dr du / r^2
  if (!doForces) return;
  du /= r2;
  for (int k=0; k<3; k++) {
    fi[k] += dr[k]*du;
    fj[k] -= dr[k]*du;
  }
}

void PotentialMasterCell::handleComputeOne(Potential* pij, const double *ri, const double *rj, const int jAtom, double& uTot, double rc2) {
  double r2 = 0;
  for (int k=0; k<3; k++) {
    double dr = ri[k]-rj[k];
    r2 += dr*dr;
  }
  if (r2 > rc2) return;
  double uij = pij->u(r2);
  uAtomsChanged.push_back(jAtom);
  duAtom[0] += 0.5*uij;
  duAtom.push_back(0.5*uij);
  uTot += uij;
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
  double rjp[3];
  for (int i=0; i<numAtoms; i++) {
    uAtom[i] = 0;
    if (doForces) for (int k=0; k<3; k++) force[i][k] = 0;
  }
  for (int iAtom=0; iAtom<numAtoms; iAtom++) {
    int iMolecule = iAtom;
    vector<int> *iBondedAtoms = nullptr;
    if (!pureAtoms) {
      iMolecule = box.getMolecule(iAtom);
      if (!rigidMolecules) {
        int iSpecies, iLastAtom, iFirstAtom;
        box.getMoleculeInfo(iMolecule, iSpecies, iFirstAtom, iLastAtom);
        int iChildIndex = iAtom-iFirstAtom;
        iBondedAtoms = &bondedAtoms[iSpecies][iChildIndex];
      }
    }
    const int iType = box.getAtomType(iAtom);
    const double *iCutoffs = pairCutoffs[iType];
    Potential** iPotentials = pairPotentials[iType];
    const double *ri = box.getAtomPosition(iAtom);
    double *fi = doForces ? force[iAtom] : nullptr;
    int jAtom=iAtom;
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      if (checkSkip(jAtom, iMolecule, iBondedAtoms)) continue;
      const double *rj = box.getAtomPosition(jAtom);
      const int jType = box.getAtomType(jAtom);
      handleComputeAll(iAtom, jAtom, ri, rj, iPotentials[jType], uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, iCutoffs[jType], doForces);
    }
    const int iCell = atomCell[iAtom];
    for (vector<int>::iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      double *jbo = boxOffsets[jCell];
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        if (checkSkip(jAtom, iMolecule, iBondedAtoms)) continue;
        const int jType = box.getAtomType(jAtom);
        const double *rj = box.getAtomPosition(jAtom);
        for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
        handleComputeAll(iAtom, jAtom, ri, rjp, iPotentials[jType], uAtom[iAtom], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, iCutoffs[jType], doForces);
      }
    }
  }
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}

void PotentialMasterCell::computeOne(const int iAtom, const double *ri, double &u1, const bool isTrial) {
  u1 = 0;
  const int iType = box.getAtomType(iAtom);
  const double *iCutoffs = pairCutoffs[iType];
  Potential** iPotentials = pairPotentials[iType];
  double rjp[3];

  const int iCell = isTrial ? cellForCoord(ri) : atomCell[iAtom];

  uAtomsChanged.resize(1);
  duAtom.resize(1);
  uAtomsChanged[0] = iAtom;
  duAtom[0] = 0;
  int iMolecule = iAtom;
  vector<int> *iBondedAtoms = nullptr;
  if (!pureAtoms) {
    iMolecule = box.getMolecule(iAtom);
    if (!rigidMolecules) {
      int iSpecies, iLastAtom, iFirstAtom;
      box.getMoleculeInfo(iMolecule, iSpecies, iFirstAtom, iLastAtom);
      int iChildIndex = iAtom-iFirstAtom;
      iBondedAtoms = &bondedAtoms[iSpecies][iChildIndex];
    }
  }
  for (int jAtom = cellLastAtom[iCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
    if (jAtom!=iAtom) {
      if (checkSkip(jAtom, iMolecule, iBondedAtoms)) continue;
      const int jType = box.getAtomType(jAtom);
      const double *rj = box.getAtomPosition(jAtom);
      handleComputeOne(iPotentials[jType], ri, rj, jAtom, u1, iCutoffs[jType]);
    }
  }

  for (vector<int>::iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
    int jCell = iCell + *it;
    double *jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      if (checkSkip(jAtom, iMolecule, iBondedAtoms)) continue;
      const double *rj = box.getAtomPosition(jAtom);
      for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
      const int jType = box.getAtomType(jAtom);
      handleComputeOne(iPotentials[jType], ri, rjp, jAtom, u1, iCutoffs[jType]);
    }
    jCell = iCell - *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      if (checkSkip(jAtom, iMolecule, iBondedAtoms)) continue;
      const double *rj = box.getAtomPosition(jAtom);
      for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
      const int jType = box.getAtomType(jAtom);
      handleComputeOne(iPotentials[jType], ri, rjp, jAtom, u1, iCutoffs[jType]);
    }
  }
}

void PotentialMasterCell::removeAtom(int iAtom) {
  PotentialMaster::removeAtom(iAtom);
  int oldCell = atomCell[iAtom];
  if (oldCell>-1) {
    // delete from old cell
    int j = cellLastAtom[oldCell];
    if (j == iAtom) {
      cellLastAtom[oldCell] = cellNextAtom[iAtom];
    }
    else {
      while (cellNextAtom[j] != iAtom) {
        j = cellNextAtom[j];
      }
      cellNextAtom[j] = cellNextAtom[iAtom];
    }
  }
  moveAtomIndex(box.getNumAtoms()-1, iAtom);
  numAtoms--;
}

void PotentialMasterCell::moveAtomIndex(int oldIndex, int newIndex) {
  if (oldIndex==newIndex) return;
  int cell = atomCell[oldIndex];
  //printf("%d was in %d\n", oldIndex, cell);
  atomCell[newIndex] = cell;
  cellNextAtom[newIndex] = cellNextAtom[oldIndex];
  //printf("%d now points to %d\n", newIndex, cellNextAtom[newIndex]);
  if (cell==-1) return;
  int j = cellLastAtom[cell];
  if (j == oldIndex) {
    cellLastAtom[cell] = newIndex;
    return;
  }
  while (cellNextAtom[j] != oldIndex) j = cellNextAtom[j];
  cellNextAtom[j] = newIndex;
}

void PotentialMasterCell::newAtom() {
  // don't actually increment numAtoms... updateAtom will do that
  PotentialMaster::newAtom();
  updateAtom(numAtoms);
  for (int i=0; i<numAtoms; i++) {
    if (atomCell[i] < 0) {
      printf("oops %d in %d after newAtom\n", i, atomCell[i]);
    }
  }
}
