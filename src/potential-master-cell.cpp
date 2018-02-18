#include <stdio.h>
#include <iostream>
#include "alloc2d.h"
#include "potential-master.h"

PotentialMasterCell::PotentialMasterCell(Potential& p2, Box& box, double pRange, int cRange) : PotentialMaster(p2, box), range(pRange), cellRange(cRange), rawBoxOffsets(nullptr), boxOffsets(nullptr) {
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
#define mapValue (dx+(dy+dz*numCells[1])*numCells[0])
#define checkD if (dz==0) { \
                if (dy<1) { \
                 if (dx<1) continue; \
                } \
                else if (dx<0) { \
                 continue; \
                } \
               } \

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

void PotentialMasterCell::init() {
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

  int numOffsets = 0;
  switch (cellRange) {
    case 1:
      numOffsets = 13;
      break;
    case 2:
      numOffsets = 62;
      break;
    case 3:
      numOffsets = 155;
      break;
  }
  cellOffsets.resize(numOffsets);
  int o[7] = {0, 1, -1, 2, -2, 3, -3};
  int oz[4] = {0, 1, 2, 3};
  int dCell = 0;
  for (int iz=0; iz<2; iz++) {
    for (int iy=0; iy<3; iy++) {
      for (int ix=0; ix<3; ix++) {
        checkD
        cellOffsets[dCell] = mapValue;
        dCell++;
      }
    }
  }
  if (cellRange>1) {
    for (int iz=0; iz<3; iz++) {
      for (int iy=0; iy<5; iy++) {
        int ixMin = (iy<3 && iz!=2) ? 3 : 0;
        for (int ix=ixMin; ix<5; ix++) {
          checkD
          cellOffsets[dCell] = mapValue;
          dCell++;
        }
      }
    }
    if (cellRange>2) {
      for (int iz=0; iz<4; iz++) {
        for (int iy=0; iy<7; iy++) {
          if (iy>5 && iz==4) {
            for (int ix=0; ix<3; ix++) {
              cellOffsets[dCell] = mapValue;
              dCell++;
            }
            continue;
          }
          if (abs(dy)+dz==5) {  // 3+2
            for (int ix=0; ix<5; ix++) {
              cellOffsets[dCell] = mapValue;
              dCell++;
            }
            continue;
          }
          int ixMin = (iy<5 && iz!=3) ? 5 : 0;
          for (int ix=ixMin; ixMin<7; ixMin++) {
            checkD;
            cellOffsets[dCell] = mapValue;
            dCell++;
          }
        }
      }
    }
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

int PotentialMasterCell::cellForCoord(double *r) {
  int cellNum = 0;
  double* bs = box.getBoxSize();
  for (int i=0; i<3; i++) {
    double x = (r[i] + boxHalf[i])/bs[i];
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
  if (range == 0 || cellRange == 0) {
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

void PotentialMasterCell::handleComputeAll(int iAtom, int jAtom, double *ri, double *rj, double &ui, double &uj, double* fi, double* fj, double& uTot, double& virialTot, double rc2, bool doForces) {
  double r2 = 0;
  double dr[3];
  for (int k=0; k<3; k++) {
    dr[k] = rj[k]-ri[k];
    r2 += dr[k]*dr[k];
  }
  if (r2 > rc2) return;
  double u, du, d2u;
  potential.u012(r2, u, du, d2u);
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

void PotentialMasterCell::handleComputeOne(double *ri, double *rj, int jAtom, double& uTot, double rc2) {
  double r2 = 0;
  for (int k=0; k<3; k++) {
    double dr = ri[k]-rj[k];
    r2 += dr*dr;
  }
  if (r2 > rc2) return;
  double uij = potential.u(r2);
  uAtomsChanged[numAtomsChanged] = jAtom;
  duAtom[0] += 0.5*uij;
  duAtom[numAtomsChanged] = 0.5*uij;
  numAtomsChanged++;
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
  int numAtoms = box.getNumAtoms();
  double uTot=0, virialTot=0;
  double rc = potential.getCutoff();
  double rc2 = rc*rc;
  double rjp[3];
  for (int i=0; i<numAtoms; i++) {
    uAtom[i] = 0;
    if (doForces) for (int k=0; k<3; k++) force[i][k] = 0;
  }
  for (int i=0; i<numAtoms; i++) {
    double *ri = box.getAtomPosition(i);
    double *fi = doForces ? force[i] : nullptr;
    int jAtom=i;
    while ((jAtom = cellNextAtom[jAtom]) > -1) {
      double *rj = box.getAtomPosition(jAtom);
      handleComputeAll(i, jAtom, ri, rj, uAtom[i], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, rc2, doForces);
    }
    int iCell = atomCell[i];
    for (vector<int>::iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
      int jCell = iCell + *it;
      double *jbo = boxOffsets[jCell];
      jCell = wrapMap[jCell];
      for (jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
        double *rj = box.getAtomPosition(jAtom);
        for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
        handleComputeAll(i, jAtom, ri, rjp, uAtom[i], uAtom[jAtom], fi, doForces?force[jAtom]:nullptr, uTot, virialTot, rc2, doForces);
      }
    }
  }
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(uTot, virialTot, force);
  }
}

void PotentialMasterCell::computeOne(int iAtom, double *ri, double &u1, bool isTrial) {
  u1 = 0;
  double rc = potential.getCutoff();
  double rc2 = rc*rc;
  double rjp[3];

  int iCell = isTrial ? cellForCoord(ri) : atomCell[iAtom];
  numAtomsChanged = 1;
  uAtomsChanged[0] = iAtom;
  duAtom[0] = 0;
  for (int jAtom = cellLastAtom[iCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
    if (jAtom!=iAtom) {
      double *rj = box.getAtomPosition(jAtom);
      handleComputeOne(ri, rj, jAtom, u1, rc2);
    }
  }

  for (vector<int>::iterator it = cellOffsets.begin(); it!=cellOffsets.end(); ++it) {
    int jCell = iCell + *it;
    double *jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      double *rj = box.getAtomPosition(jAtom);
      for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
      handleComputeOne(ri, rjp, jAtom, u1, rc2);
    }
    jCell = iCell - *it;
    jbo = boxOffsets[jCell];
    jCell = wrapMap[jCell];
    for (int jAtom = cellLastAtom[jCell]; jAtom>-1; jAtom = cellNextAtom[jAtom]) {
      double *rj = box.getAtomPosition(jAtom);
      for (int k=0; k<3; k++) rjp[k] = rj[k] + jbo[k];
      handleComputeOne(ri, rjp, jAtom, u1, rc2);
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
