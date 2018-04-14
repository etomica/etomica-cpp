#pragma once

#include <vector>
#include "box.h"
#include "atom-info.h"

using namespace std;

class Box;

class Species {
  protected:
    int numAtoms, numAtomTypes;
    int* atomTypes;
    double** positions;
    double com[3];
    AtomInfo* atomInfo;
    void setup(int numAtoms, int numAtomTypes);

  public:
    Species(int numAtoms, int numAtomTypes);
    virtual ~Species();
    // to be called by SpeciesList when Species is added
    virtual void init(AtomInfo& atomInfo);
    int* getAtomTypes();
    int getNumAtoms();
    double* getAtomPosition(int iAtom);
    double* getMoleculeCOM(Box& box, int iFirstAtom, int iLastAtom);
};

class SpeciesSimple : public Species {
  private:
    double mass;

  public:
    SpeciesSimple(int numAtoms, double mass);
    virtual ~SpeciesSimple() {}
    virtual void init(AtomInfo& atomInfo);
    void setAtomPosition(int iAtom, double* iPosition);
};

class SpeciesFile : public Species {
  private:
    vector<double> typeMass;
    vector<int> types;
    vector<char*> typeSymbols;
    int typeOffset;

  public:
    SpeciesFile(const char *filename);
    virtual ~SpeciesFile();
    virtual void init(AtomInfo& atomInfo);
    void setAtomPosition(int iAtom, double* iPosition);
    int getTypeForSymbol(const char* symbol);
};

class SpeciesList {
  private:
    int nSpecies;
    Species** allSpecies;
    AtomInfo atomInfo;

  public:
    SpeciesList();
    ~SpeciesList();
    int size() const;
    int add(Species* species);
    Species* get(int iSpecies) const;
    AtomInfo& getAtomInfo();
    int getNumAtomTypes() const;
    bool isPurelyAtomic() const;
};
