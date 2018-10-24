/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class Matrix {
  protected:
    double* tmpVec;
  public:
    const int nRows, nCols;
    double **matrix;
    Matrix(int rows, int cols);
    ~Matrix();
    void setRows(double** rows);
    void transpose();
    void transform(double* vec);
    void TE(Matrix& m);
    void invert();
};
