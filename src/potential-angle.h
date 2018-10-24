/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class PotentialAngle {
  public:
    PotentialAngle() {}
    virtual ~PotentialAngle() {}
    virtual double u(double costheta) = 0;
    virtual double du(double costheta) = 0;
    virtual double d2u(double costheta) = 0;
    virtual void u012(double costheta, double &u, double &du, double &d2u) = 0;
};

class PotentialAngleHarmonic : public PotentialAngle {
  private:
    double k, theta0;

  public:
    PotentialAngleHarmonic(double k, double theta0);
    ~PotentialAngleHarmonic();
    double u(double costheta);
    double du(double costheta);
    double d2u(double costheta);
    void u012(double costheta, double &u, double &du, double &d2u);
};
