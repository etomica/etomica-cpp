#pragma once

#include "box.h"
#include "random.h"

class RotationMatrix {
  private:
    double matrix[3][3];

  public:
    RotationMatrix() {}
    ~RotationMatrix() {}
    void setSimpleAxisAngle(int iAxis, double theta);
    void invert();
    void transform(double* vec);
    void transformAbout(double* vec, double *center, Box& box);
    void randomize(Random &random);
};
