#include <limits>
#include "potential-bond.h"

/**
 * Harmonic bond, u = k dr^2
 */

PotentialBondHarmonic::PotentialBondHarmonic(double rEq, double eps) : Potential(TRUNC_NONE, std::numeric_limits<double>::infinity()), r0(rEq), epsilon(eps) {}

PotentialBondHarmonic::~PotentialBondHarmonic() {}

double PotentialBondHarmonic::ur(double r) {
  double dr = r - r0;
  return epsilon*dr*dr;
}

double PotentialBondHarmonic::u(double r2) {
  double r = sqrt(r2);
  double dr = r - r0;
  return epsilon*dr*dr;
}

double PotentialBondHarmonic::du(double r2) {
  double r = sqrt(r2);
  double dr = r - r0;
  return 2*epsilon*dr*r;
}

double PotentialBondHarmonic::d2u(double r2) {
  return 2*epsilon*r2;
}

void PotentialBondHarmonic::u012(double r2, double &u, double &du, double &d2u){
  double r = sqrt(r2);
  double dr = r - r0;
  u = epsilon*dr*dr;
  du = 2*epsilon*dr*r;
  d2u = 2*epsilon*r2;
}
