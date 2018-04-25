#include <math.h>
#include "rotation-matrix.h"
#include "vector.h"

RotationMatrix::RotationMatrix() {
  matrix[0][0] = matrix[0][1] = matrix[0][2] = 
  matrix[1][0] = matrix[1][1] = matrix[1][2] = 
  matrix[2][0] = matrix[2][1] = matrix[2][2] = 0;
}

void RotationMatrix::setRows(double* row0, double* row1, double* row2) {
  matrix[0][0] = row0[0]; matrix[0][1] = row0[1]; matrix[0][2] = row0[2];
  matrix[1][0] = row1[0]; matrix[1][1] = row1[1]; matrix[1][2] = row1[2];
  matrix[2][0] = row2[0]; matrix[2][1] = row2[1]; matrix[2][2] = row2[2];
}

void RotationMatrix::setSimpleAxisAngle(int iAxis, double theta) {
  double st = sin(theta);
  double ct = cos(theta);
  for (int i=0; i<3; i++) {
    matrix[i][i] = i==iAxis ? 1 : ct;
    for (int j=i+1; j<3; j++) {
      if (i==iAxis || j==iAxis) {
        matrix[i][j] = matrix[j][i] = 0;
        continue;
      }
      matrix[i][j] = -st;
      matrix[j][i] = st;
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
  for (int i=0; i<3; i++) vec[i] = newVec[i];
}

void RotationMatrix::transformAbout(double* vec, double* center, Box& box) {
  double dr[3];
  for (int i=0; i<3; i++) dr[i] = vec[i] - center[i];
  box.nearestImage(dr);
  for (int i=0; i<3; i++) {
    vec[i] = center[i];
    for (int j=0; j<3; j++) {
      vec[i] += matrix[i][j]*dr[j];
    }
  }
  box.nearestImage(vec);
}

void RotationMatrix::transpose() {
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

void RotationMatrix::TE(RotationMatrix& m) {
  for (int i=0; i<3; i++) {
    double x = matrix[i][0] * m.matrix[0][0]
             + matrix[i][1] * m.matrix[1][0]
             + matrix[i][2] * m.matrix[2][0];
    double y = matrix[i][0] * m.matrix[0][1]
             + matrix[i][1] * m.matrix[1][1]
             + matrix[i][2] * m.matrix[2][1];
    double z = matrix[i][0] * m.matrix[0][2]
             + matrix[i][1] * m.matrix[1][2]
             + matrix[i][2] * m.matrix[2][2];
    matrix[i][0] = x;
    matrix[i][1] = y;
    matrix[i][2] = z;
  }
}

double RotationMatrix::determinant() {
  return matrix[0][0]*matrix[1][1]*matrix[2][2]
        -matrix[0][0]*matrix[1][2]*matrix[2][1]
        -matrix[0][1]*matrix[1][0]*matrix[2][2]
        +matrix[0][2]*matrix[1][2]*matrix[2][1]
        +matrix[0][1]*matrix[1][2]*matrix[2][0]
        -matrix[0][2]*matrix[1][1]*matrix[2][1];
}

void RotationMatrix::invert() {
  double t00 = matrix[0][0], t01 = matrix[0][1], t02 = matrix[0][2];
  double t10 = matrix[1][0], t11 = matrix[1][1], t12 = matrix[1][2];
  double t20 = matrix[2][0], t21 = matrix[2][1], t22 = matrix[2][2];
  double det = determinant();
  matrix[0][0] = (t11*t22-t12*t21)/det; 
  matrix[0][1] = -(t01*t22-t02*t21)/det;
  matrix[0][2] = (t01*t12-t02*t11)/det;
  matrix[1][0] = -(t10*t22-t12*t20)/det; 
  matrix[1][1] = (t00*t22-t02*t20)/det;
  matrix[1][2] = -(t00*t12-t02*t10)/det;
  matrix[2][0] = (t10*t21-t11*t20)/det;
  matrix[2][1] = -(t00*t21-t01*t20)/det;
  matrix[2][2] = (t00*t11-t01*t10)/det; 
}
