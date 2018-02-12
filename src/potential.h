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
    void setCutoff(double rc);
    double getCutoff();
    void setTruncationType(int tt);
    int getTruncationType();
};

class PotentialLJ: public Potential {
  public:
    PotentialLJ(int tt, double rc);
    ~PotentialLJ() {}
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
};

class PotentialSS: public Potential {
  private:
    int exponent;
    double rpow(double r2);
  public:
    PotentialSS(int p, int tt, double rc);
    ~PotentialSS() {}
    double ur(double r);
    double u(double r2);
    double du(double r2);
    double d2u(double r2);
    void u012(double r2, double &u, double &du, double &d2u);
};

class PotentialWCA: public PotentialLJ {
  public:
    PotentialWCA();
    ~PotentialWCA() {}
};
