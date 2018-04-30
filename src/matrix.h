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
