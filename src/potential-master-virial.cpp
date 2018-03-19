#include "potential-master.h"

PotentialMasterVirial::PotentialMasterVirial(const SpeciesList &sl, Box &box) : PotentialMaster(sl,box) {}

void PotentialMasterVirial::computeAtoms(const int* iAtomList, const int nAtoms, double &energy) {
  double dr[3];
  energy = 0;
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
      for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
      double r2 = 0;
      for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
      if (r2 > iCutoffs[jType]) continue;
      double uij = iPotentials[jType]->u(r2);
      energy += uij;
    }
  }
}

void PotentialMasterVirial::computeMolecules(const int* iMoleculeList, const int nMolecules, double &energy) {
  double dr[3];
  energy = 0;
  for (int i=0; i<nMolecules-1; i++) {
    int iMolecule = iMoleculeList[i];
    int iSpecies, iFirstAtom, iLastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iFirstAtom, iLastAtom);
    for (int iAtom=iFirstAtom; iAtom<=iLastAtom; iAtom++) {
      int iType = box.getAtomType(iAtom);
      double* iCutoffs = pairCutoffs[iType];
      Potential** iPotentials = pairPotentials[iType];
      double *ri = box.getAtomPosition(iAtom);

      for (int j=i+1; j<nMolecules; j++) {
        int jMolecule = iMoleculeList[j];
        int jSpecies, jFirstAtom, jLastAtom;
        box.getMoleculeInfo(jMolecule, jSpecies, jFirstAtom, jLastAtom);
        for (int jAtom=jFirstAtom; jAtom<=jLastAtom; jAtom++) {
          //uAtomsChangedSet.insert(iAtom);
          int jType = box.getAtomType(jAtom);
          double *rj = box.getAtomPosition(jAtom);
          for (int k=0; k<3; k++) dr[k] = rj[k]-ri[k];
          double r2 = 0;
          for (int k=0; k<3; k++) r2 += dr[k]*dr[k];
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
    if ((*it)->callFinished) (*it)->allComputeFinished(0, 0, nullptr);
  }
}

double PotentialMasterVirial::uTotalFromAtoms() {
  return 0;
}
