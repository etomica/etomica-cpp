#pragma once

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
  private:
    const double epsilon;
    const int exponent;
    double rpow(double r2);
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
