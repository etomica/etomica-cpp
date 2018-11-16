#include "ewald.h"
#include "species.h"

EwaldBase::EwaldBase(const SpeciesList &sl, Box& b) : speciesList(sl), box(b), rigidMolecules(true) {
  int numAtomTypes = sl.getNumAtomTypes();
  charges = new double[numAtomTypes];
  for (int i=0; i<numAtomTypes; i++) charges[i] = 0;
}

EwaldBase::~EwaldBase() {
  delete[] charges;
}

void EwaldBase::setCharge(int iType, double q) {
  charges[iType] = q;
}

void EwaldBase::setRigidMolecules(bool rm) {
  rigidMolecules = rm;
}
