/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <math.h>
#include <algorithm>
#include "matrix.h"
#include "vector.h"
#include "alloc2d.h"

Matrix::Matrix(int n, int m) : nRows(n), nCols(m) {
  matrix = (double**)malloc2D(n, m, sizeof(double));
  std::fill(matrix[0], matrix[0]+n*m, 0.0);
  tmpVec = new double[n];
}

Matrix::~Matrix() {
  free2D((void**)matrix);
  delete[] tmpVec;
}

void Matrix::E(Matrix& m) {
  for (int i=0; i<nRows; i++) {
    for (int j=0; j<nCols; j++) {
      matrix[i][j] = m.matrix[i][j];
    }
  }
}

void Matrix::setRows(double** rows) {
  for (int i=0; i<nRows; i++) {
    for (int j=0; j<nCols; j++) {
      matrix[i][j] = rows[i][j];
    }
  }
}

void Matrix::transform(double* vec) {
  for (int i=0; i<nRows; i++) {
    tmpVec[i] = 0;
    for (int j=0; j<nCols; j++) {
      tmpVec[i] += matrix[i][j]*vec[j];
    }
  }
  for (int i=0; i<nRows; i++) vec[i] = tmpVec[i];
}

void Matrix::transpose() {
  for (int i=0; i<nRows-1; i++) {
    for (int j=i+1; j<nRows; j++) {
      double t = matrix[i][j];
      matrix[i][j] = matrix[j][i];
      matrix[j][i] = t;
    }
  }
}

// assumes square matrices
void Matrix::TE(Matrix& m) {
  for (int i=0; i<nRows; i++) {
    for (int j=0; j<nRows; j++) {
      tmpVec[j] = 0;
      for (int k=0; k<nRows; k++) {
        tmpVec[j] += matrix[i][k]*m.matrix[k][j];
      }
    }
    for (int j=0; j<nRows; j++) {
      matrix[i][j] = tmpVec[j];
    }
  }
}

void Matrix::invert() {
  double** M = matrix;

  // zero L
  for (int i=0; i<nRows; i++) {
    // pivot
    int prow = i;
    double pval = fabs(M[i][i]);
    for (int j=i+1; j<nRows; j++) {
      if (fabs(M[j][i]) > pval) {
        prow = j;
        pval = fabs(M[j][i]);
      }
    }
    pval = M[prow][i];
    if (prow != i) {
      // swap and normalize
      for (int j=i; j<nCols; j++) {
        double t = M[i][j];
        M[i][j] = M[prow][j]/pval;
        M[prow][j] = t;
      }
    }
    else {
      // just normalize
      for (int j=i+1; j<nCols; j++) {
        M[i][j] /= pval;
      }
    }
    for (int j=i+1; j<nRows; j++) {
      double Mij_iMii=M[j][i];
      for (int k=i+1; k<nCols; k++) {
        M[j][k]-=M[i][k]*Mij_iMii;
      }
      M[j][i]=0.0;
    }
  }

  // zero U
  for (int i=nRows-1; i>-1; i--) {
    for (int j=0; j<i; j++) {
      // this part not strictly necessary
      for (int k=i+1; k<nRows; k++) {
        M[j][k] = 0;
      }
      // this part is actually important
      for (int k=nRows; k<nCols; k++) {
        M[j][k]-=M[i][k]*M[j][i];
      }
    }
  }
}
