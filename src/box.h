#pragma once

#include "species.h"

class Box {
  protected:
    double boxHalf[3];
    double boxSize[3];

    double ***positions, ***velocities;

    int knownNumSpecies;
    int *numAtomsBySpecies, *numMoleculesBySpecies, *maxNumMoleculesBySpecies;
    int **firstAtom, **moleculeIdx;
    int **atomTypes;

    SpeciesList &speciesList;

  public:
    Box(SpeciesList &speciesList);
    virtual ~Box();

    double* getBoxSize() {return boxSize;}
    void boxSizeUpdated();
    int getNumAtoms();
    int getTotalNumMolecules();
    int getNumMolecules(int iSpecies);
    void nearestImage(double *dr);
    void initCoordinates();
    void setBoxSize(double x, double y, double z);
    void setNumMolecules(int iSpecies, int numMolecules);
    double* getAtomPosition(int iAtom);
    double* getAtomVelocity(int iAtom);
    int getAtomType(int iAtom);
    void enableVelocities();
    void getAtomInfo(int iMolecule, int iSpecies);
    int getFirstAtom(int iSpecies, int iMoleculeInSpecies);
    void getMoleculeInfo(int iMolecule, int &iSpecies, int &firstAtom, int &lastAtom);
    int getMolecule(int iAtom);
};
