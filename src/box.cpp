/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "alloc2d.h"
#include "box.h"
#include "matrix.h"
#include "vector.h"

#define BOX_HALF_TOL 0.5000000001

Box::Box(SpeciesList &sl) : positions(nullptr), velocities(nullptr), knownNumSpecies(sl.size()), numAtomsBySpecies(nullptr), maxNumAtomsBySpecies(nullptr), speciesNumAtoms(nullptr), numMoleculesBySpecies(nullptr), maxNumMoleculesBySpecies(nullptr), firstAtom(nullptr), atomTypes(nullptr), nTransformVectors(0), transformVectorCapacity(0), transformVectors(nullptr), tV2(nullptr), speciesList(sl), rectangular(true) {
  hDirty = false;
  edgeVectors = (double**)malloc2D(3, 3, sizeof(double));
  for (int i=0; i<3; i++) {
    boxSize[i] = 0;
    for (int j=0; j<3; j++) edgeVectors[i][j] = 0;
  }

  int ss = knownNumSpecies;
  numAtomsBySpecies = new int[ss];
  maxNumAtomsBySpecies = new int[ss];
  speciesNumAtoms = new int[ss];
  numMoleculesBySpecies = new int[ss];
  maxNumMoleculesBySpecies = new int[ss];
  positions = (double***)malloc(ss*sizeof(double**));
  velocities = nullptr;
  firstAtom = (int**)malloc(ss*sizeof(int*));
  atomTypes = (int**)malloc(ss*sizeof(int*));
  for (int i=0; i<ss; i++) {
    positions[i] = nullptr;
    if (velocities) velocities[i] = nullptr;
    numMoleculesBySpecies[i] = numAtomsBySpecies[i] = maxNumMoleculesBySpecies[i] = 0;
    atomTypes[i] = firstAtom[i] = nullptr;
    speciesNumAtoms[i] = speciesList.get(i)->getNumAtoms();
  }
  periodic[0] = periodic[1] = periodic[2] = true;
  h = new Matrix(3,3);
  hInv = new Matrix(3,3);
}

Box::Box(SpeciesList &sl, bool rec) : positions(nullptr), velocities(nullptr), knownNumSpecies(sl.size()), numAtomsBySpecies(nullptr), maxNumAtomsBySpecies(nullptr), speciesNumAtoms(nullptr), numMoleculesBySpecies(nullptr), maxNumMoleculesBySpecies(nullptr), firstAtom(nullptr), atomTypes(nullptr), nTransformVectors(0), transformVectorCapacity(0), transformVectors(nullptr), tV2(nullptr), speciesList(sl), rectangular(rec) {
  edgeVectors = (double**)malloc2D(3, 3, sizeof(double));
  for (int i=0; i<3; i++) {
    boxSize[i] = 0;
    for (int j=0; j<3; j++) edgeVectors[i][j] = 0;
  }
  h = new Matrix(3,3);
  hInv = new Matrix(3,3);

  int ss = knownNumSpecies;
  numAtomsBySpecies = new int[ss];
  maxNumAtomsBySpecies = new int[ss];
  speciesNumAtoms = new int[ss];
  numMoleculesBySpecies = new int[ss];
  maxNumMoleculesBySpecies = new int[ss];
  positions = (double***)malloc(ss*sizeof(double**));
  velocities = nullptr;
  firstAtom = (int**)malloc(ss*sizeof(int*));
  atomTypes = (int**)malloc(ss*sizeof(int*));
  for (int i=0; i<ss; i++) {
    positions[i] = nullptr;
    if (velocities) velocities[i] = nullptr;
    numMoleculesBySpecies[i] = numAtomsBySpecies[i] = maxNumMoleculesBySpecies[i] = 0;
    atomTypes[i] = firstAtom[i] = nullptr;
    speciesNumAtoms[i] = speciesList.get(i)->getNumAtoms();
  }
  periodic[0] = periodic[1] = periodic[2] = true;
}

Box::~Box() {
  if (knownNumSpecies==0) return;
  delete[] numAtomsBySpecies;
  delete[] numMoleculesBySpecies;
  delete[] maxNumMoleculesBySpecies;
  delete[] maxNumAtomsBySpecies;
  delete[] speciesNumAtoms;
  for (int i=0; i<knownNumSpecies; i++) {
    if (!positions[i]) continue;
    free2D((void**)positions[i]);
    if (velocities)free2D((void**)velocities[i]);
    free(firstAtom[i]);
    free(atomTypes[i]);
  }
  free(positions);
  free(velocities);
  free(firstAtom);
  free(atomTypes);
  free2D((void**)edgeVectors);
  if (!rectangular) {
    delete h;
    delete hInv;
    free2D((void**)transformVectors);
  }
}

void Box::initCoordinates() {
  int dimLeft = 3;
  int numMolecules = getTotalNumMolecules();
  int nCellsLeft = (numMolecules+3)/4;
  int numCells[3] = {0,0,0};
  while (dimLeft > 0) {
    double smin = 1e100;
    int dmin = 0;
    double product = 1.0;
    for (int idim = 0; idim < 3; idim++) {
      if (numCells[idim] > 0)
        continue;
      if (boxSize[idim] < smin) {
        smin = boxSize[idim];
        dmin = idim;
      }
      product *= boxSize[idim];
    }
    // round off except for last dimension (then round up)
    if (dimLeft > 1) {
      numCells[dmin] = (int) round(boxSize[dmin] * pow((nCellsLeft / product), 1.0 / dimLeft));
    }
    else {
      numCells[dmin] = nCellsLeft;
    }
    if (numCells[dmin] == 0) {
      numCells[dmin] = 1;
    }
    nCellsLeft = (nCellsLeft + numCells[dmin] - 1) / numCells[dmin];
    dimLeft--;
  }

  double cellSize[3];
  for (int i=0; i<3; i++) cellSize[i] = boxSize[i] / numCells[i];
  double basisFCC[4][3];
  for (int i=0; i<4; i++) {
    for (int j=0; j<3; j++) {
      basisFCC[i][j] = i==0?0:0.5;
    }
    if (i>0) basisFCC[i][i-1] = 0.0;
  }
  int ixyz[3];
  int iMolecule = 0;
  for (ixyz[0]=0; ixyz[0]<numCells[0]; ixyz[0]++) {
    for (ixyz[1]=0; ixyz[1]<numCells[1]; ixyz[1]++) {
      for (ixyz[2]=0; ixyz[2]<numCells[2]; ixyz[2]++) {
        for (int i=0; i<4; i++) {
          if (iMolecule == numMolecules) break;
          double ri[3];
          for (int j=0; j<3; j++) {
            ri[j] = (basisFCC[i][j] + ixyz[j] - 0.5*numCells[j]) * cellSize[j];
          }
          int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
          getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
          Species* s = speciesList.get(iSpecies);
          for (int jAtom=firstAtom; jAtom<=lastAtom; jAtom++) {
            double* jPos = s->getAtomPosition(jAtom-firstAtom);
            double* rj = getAtomPosition(jAtom);
            for (int k=0; k<3; k++) rj[k] = ri[k] + jPos[k];
            centralImage(rj);
          }
          iMolecule++;
        }
      }
    }
  }
}

void Box::setNumMolecules(int iSpecies, int n) {
  if (iSpecies<0 || iSpecies>=knownNumSpecies) {
    fprintf(stderr, "unknown species id %d\n", iSpecies);
    abort();
  }
  Species* s = speciesList.get(iSpecies);
  int sna = s->getNumAtoms();
  int na = n*sna;
  if (n>maxNumMoleculesBySpecies[iSpecies]) {
    if (velocities) velocities[iSpecies] = (double**)realloc2D((void**)velocities[iSpecies], na, 3, sizeof(double));
    positions[iSpecies] = (double**)realloc2D((void**)positions[iSpecies], na, 3, sizeof(double));
    firstAtom[iSpecies] = (int*)realloc(firstAtom[iSpecies], n*sizeof(int));
    atomTypes[iSpecies] = (int*)realloc(atomTypes[iSpecies], na*sizeof(int));
    int* speciesAtomTypes = s->getAtomTypes();
    for (int i=numMoleculesBySpecies[iSpecies]; i<n; i++) {
      for (int j=i*sna; j<(i+1)*sna; j++) {
        atomTypes[iSpecies][j] = speciesAtomTypes[j-(i*sna)];
      }
    }
    maxNumMoleculesBySpecies[iSpecies] = n;
    maxNumAtomsBySpecies[iSpecies] = na;
  }
  int fa = 0;
  for (int jSpecies=0; jSpecies<knownNumSpecies; jSpecies++) {
    if (jSpecies<iSpecies) {
      fa += numAtomsBySpecies[jSpecies];
      continue;
    }
    int iStart = jSpecies==iSpecies ? numMoleculesBySpecies[iSpecies] : 0;
    for (int i=iStart; i<n; i++) {
      firstAtom[iSpecies][i] = fa + i*sna;
    }
  }
  for (int iMolecule=numMoleculesBySpecies[iSpecies]; iMolecule<n; iMolecule++) {
    // place our new molecules at origin, using nominal conformation
    for (int j=0; j<sna; j++) {
      double* jPos = s->getAtomPosition(j);
      double* rj = positions[iSpecies][sna*iMolecule + j];
      rj[0] = jPos[0]; rj[1] = jPos[1]; rj[2] = jPos[2];
    }
  }
  numMoleculesBySpecies[iSpecies] = n;
  numAtomsBySpecies[iSpecies] = na;
}

double* Box::getAtomVelocity(int i) {
  int idx = i, iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && idx >= maxNumAtomsBySpecies[iSpecies]; iSpecies++) {
    idx -= maxNumAtomsBySpecies[iSpecies];
  }
#ifdef DEBUG
  if (idx>=maxNumAtomsBySpecies[iSpecies]) {
    printf("getAtomVecloities oops i %d is more atoms than I have\n", i);
    abort();
  }
#endif
  return velocities[0][i];
}

int Box::getGlobalMoleculeIndex(int iSpecies, int iMoleculeInSpecies) {
#ifdef DEBUG
  if (iSpecies>=knownNumSpecies) {
    printf("getGlobalMoleculeIndex oops i %d is more species than I have\n", iSpecies);
    abort();
  }
#endif
  int t = 0;
  for (int jSpecies=0; jSpecies<iSpecies; jSpecies++) {
    t += numMoleculesBySpecies[jSpecies];
  }
  return t + iMoleculeInSpecies;
}

void Box::getMoleculeInfoAtom(int iAtom, int &iMoleculeInSpecies, int &iSpecies, int &iFirstAtom) {
#ifdef DEBUG
  if (iAtom>=getNumAtoms()) {
    printf("getMoleculeInfoAtom oops i %d is more atoms than I have\n", iAtom);
    abort();
  }
#endif
  iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && iAtom >= numAtomsBySpecies[iSpecies]; iSpecies++) {
    iAtom -= numAtomsBySpecies[iSpecies];
  }
  iMoleculeInSpecies = iAtom/speciesNumAtoms[iSpecies];
  iFirstAtom = firstAtom[iSpecies][iMoleculeInSpecies];
}

void Box::getMoleculeInfo(int iMolecule, int &iSpecies, int &iMoleculeInSpecies, int &fa, int &la) {
  iMoleculeInSpecies = iMolecule;
  iSpecies = 0;
  for ( ; iSpecies<knownNumSpecies-1 && iMoleculeInSpecies >= numMoleculesBySpecies[iSpecies]; iSpecies++) {
    iMoleculeInSpecies -= numMoleculesBySpecies[iSpecies];
  }
#ifdef DEBUG
  if (iMoleculeInSpecies>=numMoleculesBySpecies[iSpecies]) {
    printf("gFA oops i %d is more molecules than I have\n", iMolecule);
    abort();
  }
#endif
  fa = firstAtom[iSpecies][iMoleculeInSpecies];
  la = fa + speciesNumAtoms[iSpecies] - 1;
}

void Box::testTransformVector(double* t) {
  for (int i=0; i<nTransformVectors; i++) {
    double dot = 0.5*fabs(Vector::dot(t, transformVectors[i]))/tV2[i];
    if (dot > 1-BOX_HALF_TOL) return;
  }
  nTransformVectors++;
  if (transformVectorCapacity < nTransformVectors) {
    transformVectors = (double**)realloc2D((void**)transformVectors, nTransformVectors, 3, sizeof(double));
    tV2 = (double*)realloc(tV2, nTransformVectors*sizeof(double));
    transformVectorCapacity = nTransformVectors;
  }
  double sq = 0;
  for (int k=0; k<3; k++) {
    transformVectors[nTransformVectors-1][k] = t[k];
    sq += t[k]*t[k];
  }
  tV2[nTransformVectors-1] = sqrt(sq);
}

void Box::boxSizeUpdated() {
  for (int i=0; i<3; i++) boxHalf[i] = 0.5*boxSize[i];
  if (!rectangular) {
    if (volume()==0) return;
    h->setRows(edgeVectors);
    h->transpose();
    hInv->E(*h);
    hInv->invert();
    // find transform vectors for nearestImage
    if (transformVectors==nullptr) {
      transformVectors = (double**)malloc2D(3, 3, sizeof(double));
      tV2 = (double*)malloc(3*sizeof(double));
      nTransformVectors = transformVectorCapacity = 3;
    }
    for (int i=0; i<3; i++) {
      double sq = 0;
      for (int j=0; j<3; j++) {
        transformVectors[i][j] = edgeVectors[i][j];
        sq += edgeVectors[i][j]*edgeVectors[i][j];
      }
      tV2[i] = sq;
    }
    // test edge pairs
    double temp1[3];
    for (int i=0; i<3-1; i++) {
      for (int j=i+1; j<3; j++) {
        double dot = Vector::dot(edgeVectors[i], edgeVectors[j]);
        if (fabs(dot) <= 1e-10) continue;
        double s = dot < 0 ? -1 : +1;
        for (int k=0; k<3; k++) temp1[k] = edgeVectors[i][k] + s*edgeVectors[j][k];
        testTransformVector(temp1);
      }
    }

    //test edge triplets
    double dot01 = Vector::dot(edgeVectors[0], edgeVectors[1]);
    double dot02 = Vector::dot(edgeVectors[0], edgeVectors[2]);
    double dot12 = Vector::dot(edgeVectors[1], edgeVectors[2]);
    
    double sum;
    sum = dot01 + dot02 + dot12;
    if (sum < 1e-10) {
      for (int k=0; k<3; k++) temp1[k] = edgeVectors[0][k] + edgeVectors[1][k] + edgeVectors[2][k];
      testTransformVector(temp1);
    }

    sum = dot01 - dot02 - dot12;
    if (sum < 1e-10) {
      for (int k=0; k<3; k++) temp1[k] = edgeVectors[0][k] + edgeVectors[1][k] - edgeVectors[2][k];
      testTransformVector(temp1);
    }

    sum = -dot01 + dot02 - dot12; 
    if (sum < 1e-10) {
      for (int k=0; k<3; k++) temp1[k] = edgeVectors[0][k] - edgeVectors[1][k] + edgeVectors[2][k];
      testTransformVector(temp1);
    }

    sum = -dot01 - dot02 + dot12; 
    if (sum < 1e-10) {
      for (int k=0; k<3; k++) temp1[k] = edgeVectors[0][k] - edgeVectors[1][k] - edgeVectors[2][k];
      testTransformVector(temp1);
    }
  }
  else {
    hDirty = true;
  }
}

double Box::volume() {
  if (rectangular) return boxSize[0]*boxSize[1]*boxSize[2];
  double temp[3];
  Vector::cross(edgeVectors[0], edgeVectors[1], temp);
  return fabs(Vector::dot(temp, edgeVectors[2]));
}

void Box::centralImage(double *dr) {
  if (rectangular) {
    nearestImage(dr);
    return;
  }
  hInv->transform(dr);
  for (int i=0; i<3; i++) {
    if (!periodic[i]) continue;
    while (dr[i] > 0.5) dr[i]--;
    while (dr[i] < -0.5) dr[i]++;
  }
  h->transform(dr);
}

void Box::nearestImage(double *dr) {
  if (rectangular) {
    for (int i=0; i<3; i++) {
      if (!periodic[i]) continue;
      while (dr[i] > boxHalf[i]) dr[i] -= boxSize[i];
      while (dr[i] < -boxHalf[i]) dr[i] += boxSize[i];
    }
  }
  else {
    // To get out of the loop, we need to check n consecutive transformVectors
    // without applying any of them.  If we reach the end, then we wrap back around.
    for (int noTransformCount=0, i=0; noTransformCount<nTransformVectors; i++) {
      double dot = Vector::dot(transformVectors[i],dr)/tV2[i];
      if (fabs(dot) > BOX_HALF_TOL) {
        dot = round(dot);
        for (int k=0; k<3; k++) dr[k] -= dot*transformVectors[i][k];
        // We transformed, but this now counts as a non-transfomrm -- we don't need to check it again.
        noTransformCount = 1;
      }
      else {
        noTransformCount++;
      }
      if (i==nTransformVectors-1) {
        // we finished a pass, but transformed a vector.
        // make another pass (we'll stop at the vector we transformed)
        i=-1;
      }
    }
  }
}

const bool* Box::getPeriodic() {
  return periodic;
}

void Box::setPeriodic(bool x, bool y, bool z) {
  periodic[0] = x;
  periodic[1] = y;
  periodic[2] = z;
}

void Box::setBoxSize(double x, double y, double z) {
  boxSize[0] = edgeVectors[0][0] = x;
  boxSize[1] = edgeVectors[1][1] = y;
  boxSize[2] = edgeVectors[2][2] = z;
  boxSizeUpdated();
}

void Box::setEdgeVector(int i, double x, double y, double z) {
  if (rectangular) {
    fprintf(stderr, "Box must be rectangular\n");
    abort();
  }
  edgeVectors[i][0] = x;
  edgeVectors[i][1] = y;
  edgeVectors[i][2] = z;
  boxSizeUpdated();
}

void Box::scaleBoxToEdgeVectors(double *ex, double *ey, double *ez) {
  // set h with new edge vectors, leave hInv with the old ones
  for (int k=0; k<3; k++) {
    edgeVectors[0][k] = ex[k];
    edgeVectors[1][k] = ey[k];
    edgeVectors[2][k] = ez[k];
  }
  h->setRows(edgeVectors);
  h->transpose();

  int nm = getTotalNumMolecules();
  // first unwrap and move molecules
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* center = species->getMoleculeCOM(*this, iFirstAtom, iLastAtom);
    for (int jAtom=iFirstAtom; jAtom<=iLastAtom; jAtom++) {
      double* rj = getAtomPosition(jAtom);
      for (int k=0; k<3; k++) rj[k] -= center[k];
      nearestImage(rj);
    }
    hInv->transform(center);
    h->transform(center);
    for (int jAtom=iFirstAtom; jAtom<=iLastAtom; jAtom++) {
      double* rj = getAtomPosition(jAtom);
      for (int k=0; k<3; k++) rj[k] += center[k];
    }
  }
  // now change box size
  setEdgeVector(0, ex[0], ex[1], ex[2]);
  setEdgeVector(1, ey[0], ey[1], ey[2]);
  setEdgeVector(2, ez[0], ez[1], ez[2]);
  int na = getNumAtoms();
  // now wrap atoms back inside box
  for (int iAtom=0; iAtom<na; iAtom++) {
    double* ri = getAtomPosition(iAtom);
    centralImage(ri);
  }
}

void Box::scaleBoxTo(double bx, double by, double bz) {

  /*for (int iAtom=0; iAtom<getNumAtoms(); iAtom++) {
    double* ri = getAtomPosition(iAtom);
    ri[0] *= bx/boxSize[0];
    ri[1] *= by/boxSize[1];
    ri[2] *= bz/boxSize[2];
  }
  setBoxSize(bx, by, bz);
  return;*/

  double s[3] = {bx/boxSize[0], by/boxSize[1], bz/boxSize[2]};
  int nm = getTotalNumMolecules();
  // first unwrap and move molecules
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* center = species->getMoleculeCOM(*this, iFirstAtom, iLastAtom);
    for (int jAtom=iFirstAtom; jAtom<=iLastAtom; jAtom++) {
      double dr[3];
      double* rj = getAtomPosition(jAtom);
      for (int k=0; k<3; k++) dr[k] = rj[k]-center[k];
      nearestImage(dr);
      for (int k=0; k<3; k++) rj[k] = center[k]*s[k] + dr[k];
    }
  }
  // now change box size
  setBoxSize(bx, by, bz);
  int na = getNumAtoms();
  // now wrap atoms back inside box
  for (int iAtom=0; iAtom<na; iAtom++) {
    double* ri = getAtomPosition(iAtom);
    centralImage(ri);
  }
}

void Box::enableVelocities() {
  if (velocities) {
    fprintf(stderr, "velocities alread enabled\n");
    return;
  }
  velocities = (double***)malloc(speciesList.size()*sizeof(double**));
  for (int i=0; i<speciesList.size(); i++) {
    int na = maxNumMoleculesBySpecies[i]*speciesList.get(i)->getNumAtoms();
    velocities[i] = (double**)malloc2D(na, 3, sizeof(double));
  }
}

Matrix* Box::getH() {
  if (rectangular && hDirty) {
    h->setRows(edgeVectors);
    h->transpose();
    hInv->E(*h);
    hDirty = false;
  }
  return h;
}

Matrix* Box::getHInv() {
  if (rectangular && hDirty) {
    h->setRows(edgeVectors);
    h->transpose();
    hInv->E(*h);
    hDirty = false;
  }
  return hInv;
}
