#pragma once

#include "potential-master.h"
#include "box.h"
#include "random.h"
#include "rotation-matrix.h"

class PotentialMaster;
class SpeciesList;

class MCMove {
  protected:
    Box& box;
    PotentialMaster& potentialMaster;
    Random& random;
    double maxStepSize;
    long numTrials, numAccepted;
    double chiSum;
    double adjustInterval;
    int lastAdjust;
    double adjustStep, minAdjustStep;

    void init();
    void adjustStepSize();

  public:
    bool verboseAdjust;
    bool tunable;
    double stepSize;

    MCMove(Box& box, PotentialMaster& potentialMaster, Random& random, double stepSize);
    virtual ~MCMove();

    virtual bool doTrial() = 0;
    virtual double getChi(double temperature) = 0;
    void setStepSize(double stepSize);
    double getStepSize();
    virtual void acceptNotify() = 0;
    virtual void rejectNotify() = 0;
    virtual double energyChange() = 0;
    double getAcceptance();
};

class MCMoveDisplacement : public MCMove {
  private:
    double rOld[3];
    double uOld, uNew;
    int iAtom;

  public:

    MCMoveDisplacement(Box& box, PotentialMaster& potentialMaster, Random& random, double stepSize);
    ~MCMoveDisplacement();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
};

class MCMoveInsertDelete : public MCMove {
  private:
    double uOld, uNew;
    int iMolecule, xMolecule;
    int firstAtom;
    bool doInsert;
    double mu;
    const int iSpecies;
    const int numAtoms;
    RotationMatrix rotMat;

  public:

    MCMoveInsertDelete(Box& box, PotentialMaster& potentialMaster, Random& random, double mu, int iSpecies);
    ~MCMoveInsertDelete();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
    bool didInsert() {return doInsert;}
};

class MCMoveMoleculeDisplacement : public MCMove {
  private:
    double deltaR[3];
    double numOldPositions;
    double **oldPositions;
    double uOld, uNew;
    int iMolecule;
    int iAtomFirst, iAtomLast;
    int iSpecies;

  public:

    MCMoveMoleculeDisplacement(Box& box, PotentialMaster& potentialMaster, Random& random, double stepSize);
    ~MCMoveMoleculeDisplacement();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
    virtual void setSpecies(int iSpecies);
};

class MCMoveMoleculeRotate : public MCMove {
  private:
    SpeciesList& speciesList;
    int numOldPositions;
    double **oldPositions;
    RotationMatrix mat;
    double uOld, uNew;
    int iMolecule;
    int iAtomFirst, iAtomLast;
    int iSpecies, mySpecies;

  public:

    MCMoveMoleculeRotate(SpeciesList& speciesList, Box& box, PotentialMaster& potentialMaster, Random& random);
    ~MCMoveMoleculeRotate();

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
    virtual void setSpecies(int iSpecies);
};

