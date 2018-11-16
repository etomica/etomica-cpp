#include "ewald.h"
#include "species.h"

EwaldBase::EwaldBase(const SpeciesList &sl, Box& b, vector<PotentialCallback*>* pcb) : speciesList(sl), box(b), pairCallbacks(pcb), rigidMolecules(true) {
  int numAtomTypes = sl.getNumAtomTypes();
  charges = new double[numAtomTypes];
  for (int i=0; i<numAtomTypes; i++) charges[i] = 0;
}

EwaldBase::~EwaldBase() {
  delete[] charges;
}

void EwaldBase::setRigidMolecules(bool rm) {
  rigidMolecules = rm;
}
