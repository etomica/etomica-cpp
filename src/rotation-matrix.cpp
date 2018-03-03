#include <math.h>
#include "rotation-matrix.h"

void RotationMatrix::setSimpleAxisAngle(int iAxis, double theta) {
  double st = sin(theta);
  double ct = cos(theta);
  for (int i=0; i<3; i++) {
    matrix[i][i] = i==iAxis ? 1 : ct;
    for (int j=i+1; j<3; j++) {
      if (i==iAxis || j==iAxis) matrix[i][j] = matrix[j][i] = 0;
      matrix[i][j] = st;
      matrix[j][i] = -st;
    }
  }
}

void RotationMatrix::transform(double* vec) {
  double newVec[3];
  for (int i=0; i<3; i++) {
    newVec[i] = 0;
    for (int j=0; j<3; j++) {
      newVec[i] += matrix[i][j]*vec[j];
    }
  }
}

void RotationMatrix::transformAbout(double* vec, double* center, Box& box) {
  double dr[3];
  for (int i=0; i<3; i++) {
    dr[i] = vec[i] - center[i];
  }
  box.nearestImage(dr);
  for (int i=0; i<3; i++) {
    vec[i] = center[i];
    for (int j=0; j<3; j++) {
      vec[i] += matrix[i][j]*dr[j];
    }
  }
  box.nearestImage(vec);
}

void RotationMatrix::invert() {
  for (int i=0; i<2; i++) {
    for (int j=i+1; j<3; j++) {
      double t = matrix[i][j];
      matrix[i][j] = matrix[j][i];
      matrix[j][i] = t;
    }
  }
}
