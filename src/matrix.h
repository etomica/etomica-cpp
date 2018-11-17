/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class Matrix {
  protected:
    double* tmpVec;
    bool arrayOwner;
  public:
    const int nRows, nCols;
    double **matrix;
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, double** array);
    ~Matrix();
    void E(Matrix& m);
    void setRows(double** rows);
    void transpose();
    void transform(double* vec);
    void TE(Matrix& m);

    /**
     * If the matrix is square, this will determine the inverse and then
     * replace the matrix with its inverse.  A => A^-1
     *
     * It the matrix is augmented (with nCols>nRows), the matrix will be inverted
     * in place with the augmented part replaced with the solution.
     * A x = b with matrix A|b, then invert will replace b with x.
     */
    void invert();
};
