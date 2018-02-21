#pragma once

class Species {
  private:
    int numAtoms;

  public:
    Species(int numAtoms);
    virtual ~Species() {}
    int getNumAtoms();
};

class SpeciesList {
  private:
    int nSpecies;
    Species** allSpecies;

  public:
    SpeciesList();
    ~SpeciesList();
    int size();
    void add(Species* species);
    Species* get(int iSpecies);
};
