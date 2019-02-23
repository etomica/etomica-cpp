/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "potential-master.h"
#include "alloc2d.h"
#include "matrix.h"

CellManager::CellManager(const SpeciesList &sl, Box& b, int cRange) : box(b), rectangular(b.isRectangular()), speciesList(sl), cellRange(cRange), range(0), rawBoxOffsets(nullptr) {
}

CellManager::~CellManager() {
  boxOffsets.clear();
  if (rawBoxOffsets) free2D((void**)rawBoxOffsets);
}

#define i_cell(ax,ay,az) ((az % numCells[2]) \
                        + (ay % numCells[1]) * numCells[2] \
                        + (ax % numCells[0]) * numCells[1] * numCells[2])

int* CellManager::getNumCells() {
  return numCells;
}

int CellManager::cellForCoord(const double *r) {
  int cellNum = 0;
  const bool* periodic = box.getPeriodic();
  if (rectangular) {
    const double* bs = box.getBoxSize();
    for (int i=0; i<3; i++) {
      const double x = (r[i] + boxHalf[i])/bs[i];
      if (periodic[i]) cellNum += ((int)(cellRange + x*(numCells[i]-2*cellRange)))*jump[i];
      else cellNum += ((int)(x*numCells[i]))*jump[i];
    }
  }
  else {
    Matrix* hInv = box.getHInv();
    double rtmp[3] = {r[0], r[1], r[2]};
    hInv->transform(rtmp);
    for (int i=0; i<3; i++) {
      const double x = rtmp[i] + 0.5;
      if (periodic[i]) cellNum += ((int)(cellRange + x*(numCells[i]-2*cellRange)))*jump[i];
      else cellNum += ((int)(x*numCells[i]))*jump[i];
    }
  }
  return cellNum;
}

void CellManager::setRange(double newRange) {
  range = newRange;
}

// returns the index of cell i inside the box, given that we have nc cells in this direction
int CellManager::wrappedIndex(int i, int nc) {
  int rv = i;
  while (rv < cellRange) rv += nc - 2*cellRange;
  while (nc - rv <= cellRange) rv -= nc - 2*cellRange;
  return rv;
}

void CellManager::init() {
  double minCellSize = range/cellRange;
  int totalCells = 1;
  const double* bs = box.getBoxSize();
  const bool* periodic = box.getPeriodic();
  double xminxyz=0, xmaxxyz=0, yminxyz=0, ymaxxyz=0, zminxyz=0, zmaxxyz=0;
  const double* edgeVectors[3] = {nullptr,nullptr,nullptr};
  int numBoxCells[3];
  if (rectangular) {
    for (int i=0; i<3; i++) {
      // with cell lists, we can accommodate rc>L/2
      // we need the box to be at least the size of a cell.
      // when rc>L/2, we'll end up doing a lattice sum
      if (bs[i] < minCellSize) {
        fprintf(stderr, "box not big enough to accomodate even 1 cell (%f < %f)\n", bs[i], minCellSize);
        exit(1);
      }
      // include cellRange of padding on each side of the box
      numBoxCells[i] = numCells[i] = ((int)floor(bs[i]/minCellSize));
      if (periodic[i]) numCells[i] += cellRange*2;
      totalCells *= numCells[i];
    }
  }
  else {
    // z first
    edgeVectors[0] = box.getEdgeVector(0);
    edgeVectors[1] = box.getEdgeVector(1);
    edgeVectors[2] = box.getEdgeVector(2);
    numBoxCells[2] = numCells[2] = (int)floor(edgeVectors[2][2]/minCellSize);
    if (periodic[2]) numCells[2] += cellRange*2;
    totalCells *= numCells[2];

    // y
    double l20 = edgeVectors[2][0]/numBoxCells[2];
    double l21 = edgeVectors[2][1]/numBoxCells[2];
    double l22 = edgeVectors[2][2]/numBoxCells[2];
    double l10 = edgeVectors[1][0];
    double l11 = edgeVectors[1][1];
    // Ly + alpha Lz
    double alpha = -(l10*l20+l11*l21)/(l20*l20+l21*l21+l22*l22);
    double maxyz = sqrt((l10+alpha*l20)*(l10+alpha*l20) + (l11+alpha*l21)*(l11+alpha*l21) + (alpha*l22)*(alpha*l22));
    numBoxCells[1] = numCells[1] = (int)floor(maxyz / minCellSize);
    if (periodic[1]) numCells[1] += cellRange*2;
    l10 /= numBoxCells[1];
    l11 /= numBoxCells[1];

    // Lx + alpha Ly + beta Lz
    // alpha = ab * b + a0
    // beta = ba * a + b0
    double l00 = edgeVectors[0][0];
    double ba = -(l10*l20 + l11*l21)/(l20*l20+l21*l21+l22*l22);
    double b0 = (l00*l20)/(l20*l20+l21*l21+l22*l22);
    double ab = -(l10*l20 + l11*l21)/(l10*l10+l11*l11+l21*l21);
    double a0 = (l00*l10)/(l10*l10+l11*l11+l21*l21);
    double beta = (ba*a0 + b0)/(1-ba*ab);
    alpha = ab*beta + a0;
    double maxxyz_x = l00+alpha*l10+beta*l20;
    double maxxyz_y = alpha*l11+beta*l21;
    double maxxyz_z = beta*l22;
    double maxxyz = sqrt(maxxyz_x*maxxyz_x + maxxyz_y*maxxyz_y + maxxyz_z*maxxyz_z);
    numBoxCells[0] = numCells[0] = (int)floor(maxxyz / minCellSize);
    if (periodic[0]) numCells[0] += cellRange*2;

    zminxyz = 0;
    zmaxxyz = edgeVectors[2][2]/numBoxCells[2];
    yminxyz = 0;
    if (edgeVectors[2][1] < 0) yminxyz += edgeVectors[2][1]/numBoxCells[2];
    ymaxxyz = edgeVectors[1][1]/numBoxCells[1];
    if (edgeVectors[2][1] > 0) yminxyz += edgeVectors[2][1]/numBoxCells[2];
    xminxyz = 0;
    if (edgeVectors[1][0] < 0) xminxyz += edgeVectors[1][0]/numBoxCells[1];
    if (edgeVectors[2][0] < 0) xminxyz += edgeVectors[2][0]/numBoxCells[2];
    xmaxxyz = edgeVectors[0][0]/numBoxCells[0];
    if (edgeVectors[1][0] > 0) xmaxxyz += edgeVectors[1][0]/numBoxCells[1];
    if (edgeVectors[2][0] > 0) xmaxxyz += edgeVectors[2][0]/numBoxCells[2];
  }
  cellLastAtom.resize(totalCells);
  wrapMap.resize(totalCells);
  boxOffsets.resize(totalCells);

  cellOffsets.resize(0);
  int dCell = 0;
  int lastCellCount = 0;

  for (int icd=1; icd<=cellRange+5; icd++) {
    // we want all cells whose squared (index) distance is not more than icd
    // but exclude any cells handle in the previous passes
    int icd2 = icd*icd;
    int iz2;
    for (int iz=0; (iz2=iz*iz)<=icd2 && iz<=cellRange; iz++) {
      int izm1 = iz==0 ? 0 : (abs(iz)-1);
      int izm1Sq = izm1*izm1;
      int iy2Max = icd2-iz2;
      int iyMax = (int)sqrt(iy2Max+0.001);
      if (iyMax>cellRange) iyMax = cellRange;
      for (int iy=-iyMax; iy<=iyMax; iy++) {
        int iym1 = iy==0 ? 0 : (abs(iy)-1);
        int iym1Sq = iym1*iym1;
        int iy2 = iy*iy;
        int ix2Max = iy2Max-iy2;
        int ixMax = (int)sqrt(ix2Max+0.001);
        if (ixMax>cellRange) ixMax = cellRange;
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
          if (rectangular) {
            if (ixm1Sq+iym1Sq+izm1Sq >= cellRange*cellRange) continue;
          }
          else {
            double dx = ix*edgeVectors[0][0]/numBoxCells[0] + iy*edgeVectors[1][0]/numBoxCells[1] + iz*edgeVectors[2][0]/numBoxCells[2];
            double min1 = xminxyz;
            double max1 = xmaxxyz;
            double min2 = xminxyz + dx;
            double max2 = xmaxxyz + dx;
            if ((max1-min2)*(max2-min1) > 0) dx = 0;
            else if (dx > 0) dx = min2 - max1;
            else dx = min1 - max2;
            double dy = iy*edgeVectors[1][1]/numBoxCells[1] + iz*edgeVectors[2][1]/numBoxCells[2];
            min1 = yminxyz;
            max1 = ymaxxyz;
            min2 = yminxyz + dy;
            max2 = ymaxxyz + dy;
            if ((max1-min2)*(max2-min1) > 0) dy = 0;
            else if (dy > 0) dy = min2 - max1;
            else dy = min1 - max2;
            double dz = iz*edgeVectors[2][2]/numBoxCells[2];
            min1 = zminxyz;
            max1 = zmaxxyz;
            min2 = zminxyz + dz;
            max2 = zmaxxyz + dz;
            if ((max1-min2)*(max2-min1) > 0) dz = 0;
            else if (dz > 0) dz = min2 - max1;
            else dz = min1 - max2;
            if (dx*dx + dy*dy + dz*dz > range*range) continue;
          }
          int mv = iz+(iy+ix*numCells[1])*numCells[2];
          cellOffsets.push_back(mv);
          dCell++;
        }
      }
    }
    if (dCell==lastCellCount) break;
    lastCellCount = dCell;
  }

  // shouldn't this be (cellRange / (numCells - 2 cellRange)) ?????
  int xboRange = periodic[0] ? (numCells[0] - cellRange - 1)/(numCells[0] - 2*cellRange) : 0;
  int yboRange = periodic[1] ? (numCells[1] - cellRange - 1)/(numCells[1] - 2*cellRange) : 0;
  int zboRange = periodic[2] ? (numCells[2] - cellRange - 1)/(numCells[2] - 2*cellRange) : 0;
  int nx = (2*xboRange+1);
  int ny = (2*yboRange+1);
  int nz = (2*zboRange+1);
  rawBoxOffsets = (double**)realloc2D((void**)rawBoxOffsets, nx*ny*nz, 3, sizeof(double));
  numRawBoxOffsets = nx*ny*nz;
  for (int ix=-xboRange; ix<=xboRange; ix++) {
    for (int iy=-yboRange; iy<=yboRange; iy++) {
      for (int iz=-zboRange; iz<=zboRange; iz++) {
        int idx = (ix+xboRange)*ny*nz+(iy+yboRange)*nz+(iz+zboRange);
        if (rectangular) {
          rawBoxOffsets[idx][0] = ix*bs[0];
          rawBoxOffsets[idx][1] = iy*bs[1];
          rawBoxOffsets[idx][2] = iz*bs[2];
        }
        else {
          for (int k=0; k<3; k++) rawBoxOffsets[idx][k] = ix*edgeVectors[0][k];
          for (int k=0; k<3; k++) rawBoxOffsets[idx][k] += iy*edgeVectors[1][k];
          for (int k=0; k<3; k++) rawBoxOffsets[idx][k] += iz*edgeVectors[2][k];
        }
      }
    }
  }

  for (int ix=0; ix<numCells[0]; ix++) {
    int x2 = periodic[0] ? wrappedIndex(ix, numCells[0]) : ix;
    int xbo = periodic[0] ? (ix-x2)/(numCells[0]-2*cellRange) : 0;
    for (int iy=0; iy<numCells[1]; iy++) {
      int y2 = periodic[1] ? wrappedIndex(iy, numCells[1]) : iy;
      int ybo = periodic[1] ? (iy-y2)/(numCells[1]-2*cellRange) : 0;
      for (int iz=0; iz<numCells[2]; iz++) {
        int z2 = periodic[2] ? wrappedIndex(iz, numCells[2]) : iz;
        int zbo = periodic[2] ? (iz-z2)/(numCells[2]-2*cellRange) : 0;
        int iMap = i_cell(ix, iy, iz);
        int dCell = i_cell(x2, y2, z2);
        wrapMap[iMap] = dCell;
        boxOffsets[iMap] = rawBoxOffsets[(xbo+xboRange)*ny*nz+(ybo+yboRange)*nz+(zbo+zboRange)];
      }
    }
  }

  assignCells();
}

void CellManager::assignCells() {
  const double *bs = box.getBoxSize();
  if (bs[0]*bs[1]*bs[2] == 0) {
    fprintf(stderr, "box has 0 volume, can't assign cells\n");
    return;
  }
  if (range == 0 || cellRange == 0) {
    fprintf(stderr, "range (%f) and cell range (%d) need to be non-zero\n", range, cellRange);
    return;
  }

  const int numAtoms = box.getNumAtoms();
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
    if (rectangular) {
      for (int i=0; i<3; i++) {
        double x = (r[i] + boxHalf[i])/bs[i];
        int y = ((int)(cellRange + x*(numCells[i]-2*cellRange)));
        if (y==numCells[i]-cellRange) y--;
        else if (y==cellRange-1) y++;
        cellNum += y*jump[i];
      }
    }
    else {
      Matrix* hInv = box.getHInv();
      double rtmp[3] = {r[0], r[1], r[2]};
      hInv->transform(rtmp);
      for (int i=0; i<3; i++) {
        const double x = rtmp[i] + 0.5;
        int y = ((int)(cellRange + x*(numCells[i]-2*cellRange)));
        if (y==numCells[i]-cellRange) y--;
        else if (y==cellRange-1) y++;
        cellNum += y*jump[i];
      }
    }
#ifdef DEBUG
    if (wrapMap[cellNum] != cellNum) {
      fprintf(stderr, "iCell %d wraps to %d for atom %d\n", cellNum, wrapMap[cellNum], iAtom);
      abort();
    }
#endif
    atomCell[iAtom] = cellNum;
    cellNextAtom[iAtom] = cellLastAtom[cellNum];
    cellLastAtom[cellNum] = iAtom;
  }
}

void CellManager::removeAtom(int iAtom) {
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
}

void CellManager::updateAtom(int iAtom) {
  int cellNum = 0;
  const double *bs = box.getBoxSize();
  double *r = box.getAtomPosition(iAtom);
  if (rectangular) {
    for (int i=0; i<3; i++) {
      double x = (r[i] + boxHalf[i])/bs[i];
      int y = ((int)(cellRange + x*(numCells[i]-2*cellRange)));
      if (y==numCells[i]-cellRange) y--;
      else if (y==cellRange-1) y++;
      cellNum += y*jump[i];
    }
  }
  else {
    Matrix* hInv = box.getHInv();
    double rtmp[3] = {r[0], r[1], r[2]};
    hInv->transform(rtmp);
    for (int i=0; i<3; i++) {
      const double x = rtmp[i] + 0.5;
      int y = ((int)(cellRange + x*(numCells[i]-2*cellRange)));
      if (y==numCells[i]-cellRange) y--;
      else if (y==cellRange-1) y++;
      cellNum += y*jump[i];
    }
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

void CellManager::removeMolecule(int iSpecies, int iMolecule) {
  // we need to fix cell lists by filling in the hole left by iMolecule with 
  // the data from the last molecule of iSpecies, and then shifts all other
  // species down
  int firstAtom = box.getFirstAtom(iSpecies, iMolecule);
  int speciesAtoms = speciesList.get(iSpecies)->getNumAtoms();
  int jMolecule = box.getNumMolecules(iSpecies)-1;
  int jFirstAtom = box.getFirstAtom(iSpecies, jMolecule);
  // we have a hole where our deleted atoms used to be.  fill it in with the
  // last atom from this species
  for (int j=0; j<speciesAtoms; j++) {
    // delete from cell lists
    removeAtom(firstAtom+j);
    // and replace
    moveAtomIndex(jFirstAtom+j, firstAtom+j);
  }
  // now shift for all species>iSpecies
  int numAtoms = box.getNumAtoms();
  for (int jAtom=jFirstAtom+speciesAtoms; jAtom<numAtoms; jAtom++) {
    moveAtomIndex(jAtom, jAtom-speciesAtoms);
  }
}

void CellManager::moveAtomIndex(int oldIndex, int newIndex) {
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

void CellManager::newMolecule(int iSpecies) {
  // we need to makes room for our atoms in cell lists and then shift
  // later species to make room.
  int iMolecule = box.getNumMolecules(iSpecies)-1;
  int firstAtom = box.getFirstAtom(iSpecies, iMolecule);
  int speciesAtoms = speciesList.get(iSpecies)->getNumAtoms();
  int lastAtom = firstAtom + speciesAtoms - 1;
  int numAtoms = box.getNumAtoms();
  cellNextAtom.resize(numAtoms);
  atomCell.resize(numAtoms);
  // first we have to shift uAtoms for all species>iSpecies
  for (int jAtom=numAtoms-1; jAtom>=firstAtom+speciesAtoms; jAtom--) {
    moveAtomIndex(jAtom, jAtom-speciesAtoms);
  }
  for (int jAtom=lastAtom; jAtom>=firstAtom; jAtom--) {
    cellNextAtom[jAtom] = -1;
    atomCell[jAtom] = -1;
    updateAtom(jAtom);
  }
#ifdef DEBUG
  for (int i=0; i<box.getNumAtoms(); i++) {
    if (atomCell[i] < 0) {
      printf("oops %d in %d after newAtom\n", i, atomCell[i]);
    }
  }
#endif
}
