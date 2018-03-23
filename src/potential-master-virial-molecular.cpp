#include "potential-master.h"
#include "potential-molecular.h"
#include "alloc2d.h"

PotentialMasterVirialMolecular::PotentialMasterVirialMolecular(const SpeciesList &sl, Box &box) : PotentialMasterVirial(sl,box) {
  moleculePairPotentials = (PotentialMolecular***)malloc2D(speciesList.size(), speciesList.size(), sizeof(PotentialMolecular*));
  for (int i=0; i<speciesList.size(); i++) {
    for (int j=0; j<speciesList.size(); j++) moleculePairPotentials[i][j] = nullptr;
  }
}

void PotentialMasterVirialMolecular::setMoleculePairPotential(int iSpecies, int jSpecies, PotentialMolecular* p) {
  moleculePairPotentials[iSpecies][jSpecies] = moleculePairPotentials[jSpecies][iSpecies] = p;
}

void PotentialMasterVirialMolecular::computeMolecules(const int* iMoleculeList, const int nMolecules, double &energy) {
  energy = 0;
  for (int i=0; i<nMolecules-1; i++) {
    int iMolecule = iMoleculeList[i];
    int iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iFirstAtom, iLastAtom);
    PotentialMolecular** iPotentials = moleculePairPotentials[iSpecies];
    for (int j=i+1; j<nMolecules; j++) {
      int jMolecule = iMoleculeList[j];
      int jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom;
      box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, jFirstAtom, jLastAtom);
      double uij = iPotentials[jSpecies]->u(box, iFirstAtom, jFirstAtom);
      energy += uij;
    }
  }
}

