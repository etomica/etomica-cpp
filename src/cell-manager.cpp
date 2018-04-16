#include "potential-master.h"
#include "alloc2d.h"

CellManager::CellManager(const SpeciesList &sl, Box& b, int cRange) : box(b), speciesList(sl), cellRange(cRange), range(0), rawBoxOffsets(nullptr) {
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
  const double* bs = box.getBoxSize();
  for (int i=0; i<3; i++) {
    const double x = (r[i] + boxHalf[i])/bs[i];
    cellNum += ((int)(cellRange + x*(numCells[i]-2*cellRange)))*jump[i];
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
  boxOffsets.resize(totalCells);

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

  int xboRange = (numCells[0] - cellRange - 1)/(numCells[0] - 2*cellRange);
  int yboRange = (numCells[1] - cellRange - 1)/(numCells[1] - 2*cellRange);
  int zboRange = (numCells[2] - cellRange - 1)/(numCells[2] - 2*cellRange);
  int nx = (2*xboRange+1);
  int ny = (2*yboRange+1);
  int nz = (2*zboRange+1);
  rawBoxOffsets = (double**)realloc2D((void**)rawBoxOffsets, nx*ny*nz, 3, sizeof(double));
  numRawBoxOffsets = nx*ny*nz;
  for (int ix=-xboRange; ix<=xboRange; ix++) {
    for (int iy=-yboRange; iy<=yboRange; iy++) {
      for (int iz=-zboRange; iz<=zboRange; iz++) {
        int idx = (ix+xboRange)*ny*nz+(iy+yboRange)*nz+(iz+zboRange);
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
    fprintf(stderr, "range and cell range need to be non-zero\n");
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
    for (int i=0; i<3; i++) {
      double x = (r[i] + boxHalf[i])/bs[i];
      int y = ((int)(cellRange + x*(numCells[i]-2*cellRange)));
      if (y==numCells[i]+2*cellRange) y--;
      cellNum += y*jump[i];
    }
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
  for (int i=0; i<3; i++) {
    double x = (r[i] + boxHalf[i])/bs[i];
    int y = ((int)(cellRange + x*(numCells[i]-2*cellRange)));
    if (y==numCells[i]+2*cellRange) y--;
    cellNum += y*jump[i];
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
