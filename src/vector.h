#pragma once

class Vector {
  public:
    static double dot(const double *v1, const double *v2) {
      return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    }
    static void cross(const double *v1, const double *v2, double *out) {
      out[0] = v1[1]*v2[2] - v1[2]*v2[1];
      out[1] = v1[2]*v2[0] - v1[0]*v2[2];
      out[2] = v1[0]*v2[1] - v1[1]*v2[0];
    }
};
