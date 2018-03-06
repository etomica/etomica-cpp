#include <math.h>
#include "rotation-matrix.h"
#include "vector.h"

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

void RotationMatrix::randomize(Random &random) {
  double *u = matrix[0], *v = matrix[1], *w = matrix[2];
  random.onSphere(u);
  double udotv;
  do {
    random.onSphere(v);
    udotv = Vector::dot(u,v);
  } while (fabs(udotv) < 0.99);
  v[0] -= udotv*u[0];
  v[1] -= udotv*u[1];
  v[2] -= udotv*u[2];
  double vNorm = sqrt(Vector::dot(v,v));
  v[0] /= vNorm;
  v[1] /= vNorm;
  v[2] /= vNorm;
  Vector::cross(u, v, w);
}
