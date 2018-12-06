/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "potential-master.h"

PotentialMasterVirial::PotentialMasterVirial(const SpeciesList &sl, Box &box) : PotentialMaster(sl,box,false) {}

void PotentialMasterVirial::computeAtoms(const int* iAtomList, const int nAtoms, double &energy) {
  for (int i=0; i<nAtoms-1; i++) {
    int iAtom = iAtomList[i];
    int iType = box.getAtomType(iAtom);
    double* iCutoffs = pairCutoffs[iType];
    Potential** iPotentials = pairPotentials[iType];
    double *ri = box.getAtomPosition(iAtom);

    for (int j=i+1; j<nAtoms; j++) {
      int jAtom = iAtomList[j];
      int jType = box.getAtomType(jAtom);
      double *rj = box.getAtomPosition(jAtom);
      double r2 = 0;
      for (int k=0; k<3; k++) {double dr = rj[k]-ri[k]; r2 += dr*dr;}
      if (r2 > iCutoffs[jType]) continue;
      double uij = iPotentials[jType]->u(r2);
      energy += uij;
    }
  }
}

void PotentialMasterVirial::computeMolecules(const int* iMoleculeList, const int nMolecules, double &energy) {
  energy = 0;
  if (pureAtoms) {
    computeAtoms(iMoleculeList, nMolecules, energy);
    return;
  }
  for (int i=0; i<nMolecules-1; i++) {
    int iMolecule = iMoleculeList[i];
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      int iType = box.getAtomType(iAtom);
      double* iCutoffs = pairCutoffs[iType];
      Potential** iPotentials = pairPotentials[iType];
      double *ri = box.getAtomPosition(iAtom);

      for (int j=i+1; j<nMolecules; j++) {
        int jMolecule = iMoleculeList[j];
        int jSpecies, iMoleculeInSpecies, jFirstAtom, jLastAtom;
        box.getMoleculeInfo(jMolecule, jSpecies, iMoleculeInSpecies, jFirstAtom, jLastAtom);
        for (int jAtom=jFirstAtom; jAtom<=jLastAtom; jAtom++) {
          //uAtomsChangedSet.insert(iAtom);
          int jType = box.getAtomType(jAtom);
          double *rj = box.getAtomPosition(jAtom);
          double r2 = 0;
          for (int k=0; k<3; k++) {double dr = rj[k]-ri[k]; r2 += dr*dr;}
          if (r2 > iCutoffs[jType]) continue;
          //uAtomsChangedSet.insert(jAtom);
          double uij = iPotentials[jType]->u(r2);
          //duAtom[jAtom] += 0.5*uij;
          //duAtom[iAtom] += 0.5*uij;
          energy += uij;
        }
      }
    }
  }
}

void PotentialMasterVirial::computeAll(vector<PotentialCallback*> &callbacks) {
  for (vector<PotentialCallback*>::iterator it = callbacks.begin(); it!=callbacks.end(); it++) {
    if ((*it)->callFinished) (*it)->allComputeFinished(0, 0, nullptr, nullptr);
  }
}

double PotentialMasterVirial::uTotalFromAtoms() {
  return 0;
}
