#include "species.h"

Species::Species(int na) : numAtoms(na) {
}

int Species::getNumAtoms() {
  return numAtoms;
}
