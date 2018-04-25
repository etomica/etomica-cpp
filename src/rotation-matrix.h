#pragma once

#include "box.h"
#include "random.h"

class Box;

class RotationMatrix {
  public:
    double matrix[3][3];
    RotationMatrix();
    ~RotationMatrix() {}
    void setSimpleAxisAngle(int iAxis, double theta);
    void setRows(double* row0, double* row1, double* row2);
    void transpose();
    void transform(double* vec);
    void transformAbout(double* vec, double *center, Box& box);
    void randomize(Random &random);
    void TE(RotationMatrix& m);
    double determinant();
    void invert();
};
