/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class Random;
class Box;

class RotationMatrix {
  public:
    double matrix[3][3];
    RotationMatrix();
    ~RotationMatrix() {}
    void setSimpleAxisAngle(int iAxis, double theta);
    void setAxisAngle(double* axis, double theta);
    void setRows(double* row0, double* row1, double* row2);
    void transpose();
    void transform(double* vec);
    void transformAbout(double* vec, double *center, Box& box);
    void randomize(Random &random);
    void TE(RotationMatrix& m);
    double determinant();
    void invert();
};
