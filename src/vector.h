/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
    static void Ev1Mv2(const double *v1, const double *v2, double *out)
    {
      out[0] = v1[0] - v2[0];
      out[1] = v1[1] - v2[1];
      out[2] = v1[2] - v2[2];
    }
    static void normalize(double *u)
    {
      double norm = 1. / sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
      u[0] *= norm;
      u[1] *= norm;
      u[2] *= norm;
    }
    static void PEa1Tv1(const double a1, const double *v1, double *out)
    {
      out[0] += a1 * v1[0];
      out[1] += a1 * v1[1];
      out[2] += a1 * v1[2];
    }
    static void PE(double *u, double *out)
    {
      out[0] += u[0];
      out[1] += u[1];
      out[2] += u[2];

    }
    static void E(double *u, double *out)
    {
      out[0] = u[0];
      out[1] = u[1];
      out[2] = u[2];
    }
    static void TE(const double a, double *out)
    {
      out[0] *= a;
      out[1] *= a;
      out[2] *= a;

    }
    static double squared(double *u)
    {
      return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    }

};
