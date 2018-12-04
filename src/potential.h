/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <stdio.h>
#include <math.h>

#define TRUNC_NONE 0
#define TRUNC_SIMPLE 1
#define TRUNC_LRC 2
#define TRUNC_SHIFT 3
#define TRUNC_FORCE_SHIFT 4

class Potential {
  protected:
    int truncType;
    double uShift, ufShift;
    double rCut;
    bool correctTruncation;
  public:
    Potential();
    Potential(int tt, double rc);
    virtual ~Potential() {}
    virtual void init();
    virtual double ur(double r) {return u(r*r);}
    virtual double u(double r2) {return 0;}
    virtual double du(double r2) {return 0;}
    virtual double d2u(double r2) {return 0;}
    virtual void u012(double r2, double &u, double &du, double &d2u);
    virtual void u012TC(double &u, double &du, double &d2u) {u=du=d2u=0;}
    void setCutoff(double rc);
    void setCorrectTruncation(bool doCorrection);
    double getCutoff();
    void setTruncationType(int tt);
    int getTruncationType();
};

class PotentialLJ: public Potential {
  protected:
    const double epsilon, sigma, sigma2;
  public:
    PotentialLJ(double epsilon, double sigma, int tt, double rc);
    ~PotentialLJ() {}
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
    virtual void u012TC(double &u, double &du, double &d2u);
};

class PotentialSS: public Potential {
  protected:
    const double epsilon;
    const int exponent;
    double epsrpow(double r2);
  public:
    PotentialSS(double epsilon, int p, int tt, double rc);
    ~PotentialSS() {}
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
    virtual void u012TC(double &u, double &du, double &d2u);
};

class PotentialSSfloat: public Potential {
  protected:
    const double epsilon;
    const double exponent;
    double epsrpow(double r2) { return epsilon*pow(r2, -0.5*exponent); }
  public:
    PotentialSSfloat(double epsilon, double p, int tt, double rc);
    ~PotentialSSfloat() {}
    virtual double ur(double r);
    virtual double u(double r2);
    virtual double du(double r2);
    virtual double d2u(double r2);
    virtual void u012(double r2, double &u, double &du, double &d2u);
    virtual void u012TC(double &u, double &du, double &d2u);
};

class PotentialSSfloatTab: public PotentialSS {
  protected:
    double** rpTab;
    const int nTab;
    const double xFac;
    double exponentFloat;
    double rpInterp(const double r2) {
      double x = r2*xFac;
      int idx = (int)x;
      x -= idx;
      double* irp = rpTab[idx];
      double z = irp[0] + x*(irp[1] + x*(irp[2] + x*irp[3]));
      /*double y = pow(r2, 0.5*exponentFloat);
      printf ("%f %f %f %e\n", sqrt(r2), z, y, (z-y)/y);*/
      return z;
    } 
  public:
    PotentialSSfloatTab(double epsilon, double p, int tt, double rc, int nTab);
    ~PotentialSSfloatTab();
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
    void u012TC(double &u, double &du, double &d2u);
};

class PotentialWCA: public PotentialLJ {
  public:
    PotentialWCA(double epsilon, double sigma);
    ~PotentialWCA() {}
    virtual void u012TC(double &u, double &du, double &d2u) {u=du=d2u=0;}
};

class PotentialHS: public Potential {
  private:
    double sigma, sigma2;
  public:
    PotentialHS(double sigma);
    ~PotentialHS() {}
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
};

class PotentialCharge : public Potential {
  private:
    const double qiqj;
    const double rCore;
  public:
    PotentialCharge(double qiqj, double rCore, double rc);
    virtual ~PotentialCharge() {}
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
};

class PotentialEwaldBare : public Potential {
  private:
    const double qiqj;
    const double alpha;
    const double twoosqrtpi;
  public:
    PotentialEwaldBare(double alpha, double qiqj, double rc);
    PotentialEwaldBare(double alpha, double qiqj, double rc, int truncType);
    virtual ~PotentialEwaldBare();
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
};

class PotentialEwald : public Potential {
  private:
    Potential& p;
    const double qiqj;
    const double alpha;
    const double twoosqrtpi;
  public:
    PotentialEwald(Potential& p2, double alpha, double qiqj, double rc);
    PotentialEwald(Potential& p2, double alpha, double qiqj, double rc, int truncType);
    virtual ~PotentialEwald();
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
    virtual void u012TC(double &u, double &du, double &d2u);
};

class PotentialEwald6 : public Potential {
  private:
    Potential& pShort;
    double Bij;
    double eta, eta2, eta6r;
  public:
    PotentialEwald6(Potential& pShort, double sigmai, double epsiloni, double sigmaj, double epsilonj, double eta, double rc);
    virtual ~PotentialEwald6();
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
    virtual void u012TC(double &u, double &du, double &d2u);
    // returns eta needed to yield desired accuracy (in terms of u(rc)/B)
    static double getEta(double rc, double uOverB);
};

