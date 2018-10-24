/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "potential-angle.h"
#include <math.h>

PotentialAngleHarmonic::PotentialAngleHarmonic(double kSpring, double t0) : k(kSpring), theta0(t0) {}

PotentialAngleHarmonic::~PotentialAngleHarmonic() {}

double PotentialAngleHarmonic::u(double costheta) {
  double dtheta = acos(costheta) - theta0;
  return k*dtheta*dtheta;
}

double PotentialAngleHarmonic::du(double costheta) {
  double theta = acos(costheta);
  return 2*k*theta;
}

double PotentialAngleHarmonic::d2u(double costheta) {
  return 2*k;
}

void PotentialAngleHarmonic::u012(double costheta, double &u, double &du, double &d2u) {
  double dtheta = acos(costheta) - theta0;
  u = k*dtheta*dtheta;
  du = 2*k*dtheta;
  d2u = 2*k;
}
