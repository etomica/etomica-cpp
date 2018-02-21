#include <cstddef>
#include <stdlib.h>
#include "species.h"

SpeciesList::SpeciesList() : nSpecies(0), allSpecies(nullptr) {}

SpeciesList::~SpeciesList() {
  if (allSpecies) free(allSpecies);
}

int SpeciesList::size() {
  return nSpecies;
}

void SpeciesList::add(Species* s) {
  allSpecies = (Species**)realloc(allSpecies, (nSpecies+1)*sizeof(Species*));
  allSpecies[nSpecies] = s;
  nSpecies++;
  s->init(atomInfo);
}

Species* SpeciesList::get(int i) {
  return allSpecies[i];
}
