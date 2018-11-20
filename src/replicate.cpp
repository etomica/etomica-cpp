/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#include "action.h"
#include "box.h"

void Replicate::go(Box& box, int replicates[3]) {
  SpeciesList& speciesList = box.getSpeciesList();
  int nRep = replicates[0]*replicates[1]*replicates[2];
  int nm = box.getTotalNumMolecules();
  // first unwrap
  for (int iMolecule=0; iMolecule<nm; iMolecule++) {
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    Species* species = speciesList.get(iSpecies);
    double* center = species->getMoleculeCOM(box, iFirstAtom, iLastAtom);
    for (int jAtom=iFirstAtom; jAtom<=iLastAtom; jAtom++) {
      double dr[3];
      double* rj = box.getAtomPosition(jAtom);
      for (int k=0; k<3; k++) dr[k] = rj[k]-center[k];
      box.nearestImage(dr);
      for (int k=0; k<3; k++) rj[k] = center[k] + dr[k];
    }
  }
  const double* bs = box.getBoxSize();
  double newBS[3];
  for (int k=0; k<3; k++) newBS[k] = bs[k]*replicates[k];
  box.setBoxSize(newBS[0], newBS[1], newBS[2]);
  for (int iAtom=0; iAtom<box.getNumAtoms(); iAtom++) {
    double* ri = box.getAtomPosition(iAtom);
    for (int k=0; k<3; k++) ri[k] -= 0.5*(replicates[k]-1)*bs[k];
  }

  int nOffset = 0;
  for (int iSpecies=0; iSpecies<speciesList.size(); iSpecies++) {
    Species* species = speciesList.get(iSpecies);
    int N = box.getNumMolecules(iSpecies);
    int na = species->getNumAtoms();
    box.setNumMolecules(iSpecies, N*nRep);
    for (int iAtom=nOffset; iAtom<nOffset+N*na; iAtom++) {
      double* ri = box.getAtomPosition(iAtom);
      int jAtom = iAtom;
      for (int xr=0; xr<replicates[0]; xr++) {
        for (int yr=0; yr<replicates[1]; yr++) {
          for (int zr=0; zr<replicates[2]; zr++) {
            if (xr+yr+zr==0) continue;
            jAtom += N*na;
            double* rj = box.getAtomPosition(jAtom);
            int xrr[3] = {xr,yr,zr};
            for (int k=0; k<3; k++) rj[k] = ri[k] + xrr[k]*bs[k]/replicates[k];
          }
        }
      }
    }
    nOffset += N*na*nRep;
  }
  for (int iAtom=0; iAtom<box.getNumAtoms(); iAtom++) {
    double* ri = box.getAtomPosition(iAtom);
    box.nearestImage(ri);
  }
}
