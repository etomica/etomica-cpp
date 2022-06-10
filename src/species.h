/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <stdio.h>
#include <vector>
#include "atom-info.h"
#include "rigid-constraint.h"

using namespace std;

class Box;

class Species {
  protected:
    int numAtoms, numAtomTypes;
    int* atomTypes;
    double** positions;
    double com[3];
    AtomInfo* atomInfo;
    int axisAtoms[2][2];
    vector<RigidConstraint*> rigidConstraints;
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
    void getMoleculeOrientation(Box& box, int iFirstAtom, double* direction1, double* direction2);
    virtual vector<RigidConstraint*> getRigidConstraints();
    double getMass(int iAtom);
};

class SpeciesSimple : public Species {
  private:
    vector<double> typeMass;
    vector<double> numAtomsOfType;

  public:
    SpeciesSimple(int numAtoms, int numAtomTypes);
    virtual ~SpeciesSimple() {}
    virtual void init(AtomInfo& atomInfo);
    virtual void addAtomType(double m, int aot);
    void setAtomPosition(int iAtom, double x, double y, double z);
};

class SpeciesFile : public Species {
  private:
    vector<double> typeMass;
    vector<int> types;
    vector<char*> typeSymbols;
    int typeOffset;
    vector<RigidConstraint> rigidConstraints;

    void readAtomTypes(FILE* f, const char* filename);
    char* readAtoms(FILE* f, const char* filename, vector<double*> &tmpPositions);
    char* readOrientations(FILE* f, const char* filename);
    char* readConstraints(FILE* f, const char* filename);

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
