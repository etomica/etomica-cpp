#include <stdio.h>
#include <math.h>
#include "potential.h"

Potential::Potential(int tt, double rc) {
  truncType = tt;
  rCut = rc;
  init();
}

Potential::Potential() {
  truncType = TRUNC_SIMPLE;
  rCut = 3;
  init();
}

void Potential::init() {
  uShift = 0;
  ufShift = 0;
  if (truncType == TRUNC_NONE) {
    rCut = 1.0/0.0;
  }
  else if (truncType == TRUNC_SHIFT) {
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


PotentialLJ::PotentialLJ(int tt, double rc) : Potential(tt, rc) {
  init();
}

double PotentialLJ::ur(double r) {
  if (r>rCut) return 0;
  double s2 = 1/(r*r);
  double s6 = s2*s2*s2;
  double ulj = 4*s6*(s6 - 1);
  return ulj + uShift + r*ufShift;
}

double PotentialLJ::u(double r2) {
  if (r2>rCut*rCut) return 0;
  double s6 = 1/(r2*r2*r2);
  double ulj = 4*s6*(s6 - 1);
  double u = ulj + uShift;
  if (ufShift!=0) u += ufShift*sqrt(r2);
  return u;
}

double PotentialLJ::du(double r2) {
  double s2 = 1/r2;
  double s6 = s2*s2*s2;
  double du = -48*s6*(s6 -0.5);
  du += ufShift;
  return du;
}

double PotentialLJ::d2u(double r2) {
  double s2 = 1/r2;
  double s6 = s2*s2*s2;
  double d2u = 4*12*s6*(13*s6 - 0.5*7);
  return d2u;
}

void PotentialLJ::u012(double r2, double &u, double &du, double &d2u) {
  double s2 = 1/r2;
  double s6 = s2*s2*s2;
  u = 4*s6*(s6 - 1) + uShift;
  du = -4*12*s6*(s6 - 0.5);
  if (ufShift != 0) {
    u += sqrt(r2)*ufShift;
    du += ufShift;
  }
  d2u = 4*12*s6*(13*s6 - 0.5*7);
}

PotentialSS::PotentialSS(int p, int tt, double rc) : Potential(tt, rc) {
  exponent = p;
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
  return uss;
}

double PotentialSS::ur(double r) {
  double uss = rpow(r*r);
  return uss + uShift + r*ufShift;
}

double PotentialSS::u(double r2) {
  double u = rpow(r2) + uShift;
  if (ufShift!=0) u += uShift*sqrt(r2);
  return u;
}

double PotentialSS::du(double r2) {
  double du = -exponent * rpow(r2) + ufShift;
  return du;
}

double PotentialSS::d2u(double r2) {
  return exponent*(exponent+1)*rpow(r2);
}

void PotentialSS::u012(double r2, double &u, double &du, double &d2u) {
  u = rpow(r2);
  du = -exponent*u;
  d2u = -(exponent+1)*du;
  u += uShift;
  if (ufShift != 0) {
    u += sqrt(r2)*ufShift;
    du += ufShift;
  }
}

PotentialWCA::PotentialWCA() : PotentialLJ(TRUNC_SHIFT, 1.122462048309373) {}
