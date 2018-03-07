#include <stdio.h>
#include <math.h>
#include "potential.h"

Potential::Potential(int tt, double rc) : truncType(tt), rCut(rc) {
  init();
}

Potential::Potential() : truncType(TRUNC_SIMPLE), rCut(3) {
  init();
}

void Potential::init() {
  uShift = 0;
  ufShift = 0;
  if (truncType == TRUNC_SHIFT) {
    uShift = -ur(rCut);
  }
  else if (truncType == TRUNC_FORCE_SHIFT) {
    double rdudr = du(rCut*rCut);
    ufShift = -rdudr/rCut;
    uShift = -ur(rCut);
  }
}

void Potential::setTruncationType(int tt) {
  truncType = tt;
  init();
}

int Potential::getTruncationType() {
  return truncType;
}

void Potential::setCutoff(double rc) {
  rCut = rc;
  init();
}

double Potential::getCutoff() {
  return rCut;
}

void Potential::u012(double r2, double &u, double &du, double &d2u) {
  u = this->u(r2);
  du = this->du(r2);
  d2u = this->d2u(r2);
}


PotentialLJ::PotentialLJ(double e, double s, int tt, double rc) : Potential(tt, rc), epsilon(e), sigma(s), sigma2(s*s) {
  init();
}

double PotentialLJ::ur(double r) {
  double s2 = sigma2/(r*r);
  double s6 = s2*s2*s2;
  double ulj = 4*epsilon*s6*(s6 - 1);
  return ulj + uShift + r*ufShift;
}

double PotentialLJ::u(double r2) {
  double s2 = sigma2/r2;
  double s6 = s2*s2*s2;
  double ulj = 4*epsilon*s6*(s6 - 1);
  double u = ulj + uShift;
  if (ufShift!=0) u += ufShift*sqrt(r2);
  return u;
}

double PotentialLJ::du(double r2) {
  double s2 = 1/r2;
  double s6 = s2*s2*s2;
  double du = -48*epsilon*s6*(s6 -0.5);
  if (ufShift!=0) du += ufShift*sqrt(r2);
  return du;
}

double PotentialLJ::d2u(double r2) {
  double s2 = sigma2/r2;
  double s6 = s2*s2*s2;
  double d2u = 4*12*epsilon*s6*(13*s6 - 0.5*7);
  return d2u;
}

void PotentialLJ::u012(double r2, double &u, double &du, double &d2u) {
  double s2 = sigma2/r2;
  double s6 = s2*s2*s2;
  u = 4*epsilon*s6*(s6 - 1) + uShift;
  du = -4*12*epsilon*s6*(s6 - 0.5);
  if (ufShift != 0) {
    double x = sqrt(r2)*ufShift;
    u += x;
    du += x;
  }
  d2u = 4*12*epsilon*s6*(13*s6 - 0.5*7);
}

void PotentialLJ::u012TC(double &u0, double &u1, double &du0, double &du1, double &d2u1) {
  if (truncType == TRUNC_NONE) {
    u0=u1=du0=du1=d2u1=0;
    return;
  }
  double rc3 = rCut*rCut*rCut;
  double sc = sigma/rCut;
  double sc3 = sc*sc*sc;
  double sc6 = sc3*sc3;
  double sc12 = sc6*sc6;
  // correction due to shift and force-shift
  du0 = M_PI*rc3*rCut*ufShift;
  u0 = 4*M_PI*uShift*rc3/3 + du0;
  // correction due to truncation
  u1 = 4*M_PI*4*epsilon*(sc12/(12-3) - sc6/(6-3))*rc3;
  du1 = -4*M_PI*4*epsilon*12*(sc12/(12-3) - 0.5*sc6/(6-3))*rc3;
  d2u1 = 4*M_PI*4*epsilon*12*(13*sc12/(12-3) - 0.5*7*sc6/(6-3))*rc3;
}

PotentialSS::PotentialSS(double e, int p, int tt, double rc) : Potential(tt, rc), epsilon(e), exponent(p) {
  init();
}

double PotentialSS::rpow(double r2) {
  double s2 = 1/r2;
  double uss;
  double s6 = 0;
  if (exponent>=6) s6 = s2*s2*s2;
  switch (exponent) {
    case 4:
      uss = s2*s2;
      break;
    case 5:
      uss = s2*s2/sqrt(r2);
      break;
    case 6:
      uss = s6;
      break;
    case 8:
      uss = s6*s2;
      break;
    case 9:
      uss = s6*s2/sqrt(r2);
      break;
    case 10:
      uss = s6*s2*s2;
      break;
    case 12:
      uss = s6*s6;
      break;
    default:
      uss = pow(s2, exponent*0.5);
  }
  return epsilon*uss;
}

double PotentialSS::ur(double r) {
  double uss = rpow(r*r);
  return epsilon*uss + uShift + r*ufShift;
}

double PotentialSS::u(double r2) {
  double u = epsilon*rpow(r2) + uShift;
  if (ufShift!=0) u += uShift*sqrt(r2);
  return u;
}

double PotentialSS::du(double r2) {
  double du = -exponent * epsilon * rpow(r2);
  if (ufShift!=0) du += ufShift*sqrt(r2);
  return du;
}

double PotentialSS::d2u(double r2) {
  return exponent*(exponent+1)*epsilon*rpow(r2);
}

void PotentialSS::u012(double r2, double &u, double &du, double &d2u) {
  u = epsilon*rpow(r2);
  du = -exponent*u;
  d2u = -(exponent+1)*du;
  u += uShift;
  if (ufShift != 0) {
    double x = sqrt(r2)*ufShift;
    u += x;
    du += x;
  }
}

void PotentialSS::u012TC(double &u0, double &u1, double &du0, double &du1, double &d2u1) {
  if (truncType == TRUNC_NONE) {
    u0=u1=du0=du1=d2u1=0;
    return;
  }
  double rc3 = rCut*rCut*rCut;
  double x = epsilon*rpow(rCut);
  // correction due to shift and force-shift
  du0 = M_PI*rc3*rCut*ufShift;
  u0 = 4*M_PI*uShift*rc3/3 + du0;
  // correction due to truncation
  u1 = 4*M_PI*epsilon*x/(exponent-3);
  du1 = -u1*exponent;
  d2u1 = -du1*(exponent+1);
}

PotentialWCA::PotentialWCA(double e, double s) : PotentialLJ(e, s, TRUNC_SHIFT, 1.122462048309373*s) {}

PotentialHS::PotentialHS(double s) : Potential(TRUNC_SIMPLE, 1), sigma(s), sigma2(s*s) {
  init();
}

double PotentialHS::ur(double r) {
  return r>sigma ? 0 : INFINITY;
}

double PotentialHS::u(double r2) {
  return r2>sigma2 ? 0 : INFINITY;
}

double PotentialHS::du(double r2) {
  return 0;
}

double PotentialHS::d2u(double r2) {
  return 0;
}

void PotentialHS::u012(double r2, double &u, double &du, double &d2u) {
  u = r2>sigma2 ? 0 : INFINITY;
  du = d2u = 0;
}
