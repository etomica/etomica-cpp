/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "potential-molecular.h"

/**
 * This potential acts like a molecular potential but returns the energy
 * a pair of atoms would have if placed at the geometric centers of the molecules.
 *
 * This potential is probably only useful for hard sphere reference for
 * Mayer sampling.
 */

PotentialMolecularAtomic::PotentialMolecularAtomic(int na, Potential& p) : PotentialMolecular(), nAtoms(na), potential(p) {
}

double PotentialMolecularAtomic::u(Box& box, int iFirstAtom, int jFirstAtom) {
  double ri[3] = {0.0,0.0,0.0};
  for (int iAtom = iFirstAtom; iAtom<iFirstAtom+nAtoms; iAtom++) {
    double* rii = box.getAtomPosition(iAtom);
    ri[0] += rii[0]; ri[1] += rii[1]; ri[2] += rii[2];
  }
  ri[0] /= nAtoms; ri[1] /= nAtoms; ri[2] /= nAtoms;
  double rj[3] = {0.0,0.0,0.0};
  for (int jAtom = jFirstAtom; jAtom<jFirstAtom+nAtoms; jAtom++) {
    double* rjj = box.getAtomPosition(jAtom);
    rj[0] += rjj[0]; rj[1] += rjj[1]; rj[2] += rjj[2];
  }
  rj[0] /= nAtoms; rj[1] /= nAtoms; rj[2] /= nAtoms;
  double dr, r2=0;
  dr = ri[0]-rj[0]; r2 += dr*dr;
  dr = ri[1]-rj[1]; r2 += dr*dr;
  dr = ri[2]-rj[2]; r2 += dr*dr;
  return potential.u(r2);
}
