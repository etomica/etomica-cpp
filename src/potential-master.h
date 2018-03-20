#pragma once

#include <vector>
#include <math.h>
#include <cstddef>
#include <set>
#include <algorithm>
#include "box.h"
#include "potential.h"
#include "potential-molecular.h"

using namespace std;

class PotentialCallback {
  public:
    bool callPair;
    bool callFinished;
    bool takesForces;

    PotentialCallback();
    virtual ~PotentialCallback() {}
    virtual void reset() {}
    virtual void pairCompute(int iAtom, int jAtom, double* dr, double u, double du, double d2u) {}
    virtual void allComputeFinished(double uTot, double virialTot, double** f) {}
    virtual int getNumData() {return 0;}
    virtual double* getData() {return nullptr;}
};

class PotentialCallbackPressure : public PotentialCallback {
  protected:
    Box& box;
    double temperature;
    double data[1];
  public:
    PotentialCallbackPressure(Box& box, double temperature);
    ~PotentialCallbackPressure() {}
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};
    
class PotentialCallbackHMA : public PotentialCallback {
  protected:
    Box& box;
    double temperature;
    double Pharm;
    double data[2];
    double** latticePositions;
  public:
    PotentialCallbackHMA(Box& box, double temperature, double Pharm);
    ~PotentialCallbackHMA();
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};

class CellManager {
  public:
    Box &box;
    const SpeciesList& speciesList;
    int cellRange;
    double range;
    double boxHalf[3];
    int numCells[3];
    vector<int> cellNextAtom;
    vector<int> atomCell;
    vector<int> cellLastAtom;
    int jump[3];
    vector<int> cellOffsets;
    vector<int> wrapMap;
    double** rawBoxOffsets;
    vector<double*> boxOffsets;
    int wrappedIndex(int i, int nc);
    void moveAtomIndex(int oldIndex, int newIndex);

    CellManager(const SpeciesList &sl, Box& box, int cRange);
    ~CellManager();
    void setRange(double newRange);
    void init();
    void updateAtom(int iAtom);
    void newMolecule(int iSpecies);
    void removeAtom(int iAtom);
    void removeMolecule(int iSpecies, int iMolecule);
    void assignCells();
    int cellForCoord(const double *r);
    int* getNumCells();
};

class PotentialMaster {
  protected:
    const SpeciesList& speciesList;
    Potential*** pairPotentials;
    int* numAtomsByType;
    double** pairCutoffs;
    Box& box;
    vector<double> uAtom;
    vector<double> duAtom;
    bool duAtomSingle, duAtomMulti;
    vector<int> uAtomsChanged;
    double** force;
    vector<PotentialCallback*> pairCallbacks;
    const int numAtomTypes;
    vector<vector<int*> > *bondedPairs;
    vector<int> **bondedAtoms;
    vector<Potential*> *bondedPotentials;
    const bool pureAtoms;
    bool rigidMolecules;
    bool doTruncationCorrection, doSingleTruncationCorrection;

    void computeOneMoleculeBonds(const int iSpecies, const int iMolecule, double &u1);
    void computeAllBonds(bool doForces, double &uTot, double &virialTot);
    void computeAllTruncationCorrection(double &uTot, double &virialTot);
    double computeOneTruncationCorrection(const int iAtom);
    inline bool checkSkip(int jAtom, int iSpecies, int iMolecule, vector<int> *iBondedAtoms) {
      if (pureAtoms) return false;
      int jMolecule, jFirstAtom, jSpecies;
      box.getMoleculeInfoAtom(jAtom, jMolecule, jSpecies, jFirstAtom);
      if (rigidMolecules) return iSpecies==jSpecies && iMolecule == jMolecule;
      return binary_search(iBondedAtoms->begin(), iBondedAtoms->end(), jAtom-jFirstAtom);
    }
    virtual void computeOneInternal(const int iAtom, const double *ri, double &u1, const int iSpecies, const int iMolecule, const int iFirstAtom);

  public:
    PotentialMaster(const SpeciesList &speciesList, Box& box);
    virtual ~PotentialMaster();
    Box& getBox();
    void setDoTruncationCorrection(bool doCorrection);
    void setDoSingleTruncationCorrection(bool doCorrection);
    virtual void setPairPotential(int iType, int jType, Potential* pij);
    void setBondPotential(int iSpecies, vector<int*> &bondedPairs, Potential *pBond);
    // compute for the whole box
    virtual void computeAll(vector<PotentialCallback*> &callbacks);
    // energy of one atom with the whole box
    virtual void computeOne(const int iAtom, double &energy);
    // energy of one molecule with the whole box (including itself)
    virtual void computeOneMolecule(int iMolecule, double &energy);
    virtual void updateAtom(int iAtom) {}
    virtual void newMolecule(int iSpecies);
    virtual void removeMolecule(int iSpecies, int iMolecule);
    double oldEnergy(int iAtom);
    double oldMoleculeEnergy(int iAtom);
    void resetAtomDU();
    void processAtomU(int coeff);
    void addCallback(PotentialCallback* pcb);
    virtual double uTotalFromAtoms();
};

class PotentialMasterCell : public PotentialMaster {
  protected:
    CellManager cellManager;
    const int cellRange;
    const vector<int> &cellNextAtom;
    const vector<int> &atomCell;
    const vector<int> &cellLastAtom;
    const vector<int> &cellOffsets;
    const vector<int> &wrapMap;
    const vector<double*> &boxOffsets;
    bool lsNeeded;

    void handleComputeOne(Potential* pij, const double *ri, const double *rj, const double* jbo, const int iAtom, const int jAtom, double& uTot, double rc2) {
      double dx = ri[0]-(rj[0]+jbo[0]);
      double dy = ri[1]-(rj[1]+jbo[1]);
      double dz = ri[2]-(rj[2]+jbo[2]);
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 > rc2) return;
      double uij = pij->u(r2);
      if (duAtomSingle) {
        uAtomsChanged.push_back(jAtom);
        duAtom[0] += 0.5*uij;
        duAtom.push_back(0.5*uij);
      }
      else {
        if (duAtom[jAtom] == 0) {
          uAtomsChanged.push_back(jAtom);
        }
        duAtom[iAtom] += 0.5*uij;
        duAtom[jAtom] += 0.5*uij;
      }
      uTot += uij;
    }
    void handleComputeAll(int iAtom, int jAtom, const double *ri, const double *rj, const double *jbo, Potential* pij, double &ui, double &uj, double* fi, double* fj, double& uTot, double& virialTot, double rc2, bool doForces) {
      double dr[3];
      dr[0] = (rj[0]+jbo[0])-ri[0];
      dr[1] = (rj[1]+jbo[1])-ri[1];
      dr[2] = (rj[2]+jbo[2])-ri[2];
      double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      if (r2 > rc2) return;
      double u, du, d2u;
      pij->u012(r2, u, du, d2u);
      ui += 0.5*u;
      uj += 0.5*u;

      uTot += u;
      virialTot += du;
      for (vector<PotentialCallback*>::iterator it = pairCallbacks.begin(); it!=pairCallbacks.end(); it++) {
        (*it)->pairCompute(iAtom, jAtom, dr, u, du, d2u);
      }

      // f0 = dr du / r^2
      if (!doForces) return;
      du /= r2;
      for (int k=0; k<3; k++) {
        fi[k] += dr[k]*du;
        fj[k] -= dr[k]*du;
      }
    }

    virtual void computeOneInternal(const int iAtom, const double *ri, double &energy, const int iSpecies, const int iMolecule, const int iFirstAtom);
  public:
    PotentialMasterCell(const SpeciesList &speciesList, Box& box, int cellRange);
    ~PotentialMasterCell();
    virtual double getRange();
    virtual void setPairPotential(int iType, int jType, Potential* pij);
    // needs to be called at startup and any time the box size changes
    virtual void init();
    virtual void computeAll(vector<PotentialCallback*> &callbacks);
    virtual void updateAtom(int iAtom);
    virtual void newMolecule(int iSpecies);
    virtual void removeAtom(int iAtom);
    virtual void removeMolecule(int iSpecies, int iMolecule);
    int* getNumCells();
};

class PotentialMasterList : public PotentialMasterCell {
  protected:
    double nbrRange;
    int **nbrs;
    bool onlyUpNbrs; // standard MD only needs up.  MC or DMD needs down
    int *numAtomNbrsUp, *numAtomNbrsDn;
    int nbrsNumAtoms;
    int maxNab;
    double ***nbrBoxOffsets;
    bool forceReallocNbrs;
    double **oldAtomPositions;
    double safetyFac;
    double *maxR2, *maxR2Unsafe;

    int checkNbrPair(int iAtom, int jAtom, double *ri, double *rj, double rc2, double *jbo);
  public:
    PotentialMasterList(const SpeciesList& speciesList, Box& box, int cellRange, double nbrRange);
    ~PotentialMasterList();
    virtual double getRange();
    virtual void init();
    void reset();
    void setDoDownNbrs(bool doDown);
    void checkUpdateNbrs();
    virtual void computeAll(vector<PotentialCallback*> &callbacks);
};

class PotentialMasterVirial : public PotentialMaster {
  public:
    PotentialMasterVirial(const SpeciesList &speciesList, Box& box);
    virtual ~PotentialMasterVirial() {}
    // inter-molecular energy of a group of molecules
    void computeMolecules(const int* iMolecules, const int nMolecules, double &energy);
    void computeAtoms(const int* iAtoms, const int nAtoms, double &energy);
    void computeAll(vector<PotentialCallback*> &callbacks);
    double uTotalFromAtoms();
};

class PotentialMasterVirialMolecular : public PotentialMasterVirial {
  protected:
    PotentialMolecular*** moleculePairPotentials;

  public:
    PotentialMasterVirialMolecular(const SpeciesList &speciesList, Box& box);
    virtual ~PotentialMasterVirialMolecular() {}
    // inter-molecular energy of a group of molecules
    void setMoleculePairPotential(int iSpecies, int jSpecies, PotentialMolecular* p);
    void computeMolecules(const int* iMolecules, const int nMolecules, double &energy);
};

