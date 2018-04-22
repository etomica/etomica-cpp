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
