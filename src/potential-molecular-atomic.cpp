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

PotentialMolecularAtomic::PotentialMolecularAtomic(SpeciesList& sl, Potential& p) : PotentialMolecular(), potential(p), speciesList(sl)
{
}

double PotentialMolecularAtomic::u(Box& box, int iMolecule, int jMolecule) {
  int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;

  box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);

  double* com = speciesList.get(iSpecies)->getMoleculeCOM(box, iFirstAtom, iLastAtom);
  double ri[3] = {com[0], com[1], com[2]};
  int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;

  box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom);

  com = speciesList.get(jSpecies)->getMoleculeCOM(box, jFirstAtom, jLastAtom);

  double dr, r2=0;
  dr = ri[0]-com[0]; r2 += dr*dr;
  dr = ri[1]-com[1]; r2 += dr*dr;
  dr = ri[2]-com[2]; r2 += dr*dr;
  // printf("%f %f %f %f %f %f %f\n", ri[0], ri[1], ri[2], com[0], com[1], com[2], r2);
  return potential.u(r2);
}
