#pragma once

#include <vector>
#include <math.h>
#include "box.h"
#include "potential.h"

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
    ~PotentialCallbackHMA() {}
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};

class PotentialMaster {
  protected:
    Potential& potential;
    Box& box;
    vector<double> uAtom;
    vector<double> duAtom;
    vector<int> uAtomsChanged;
    int numAtomsChanged;
    double** force;
    vector<PotentialCallback*> pairCallbacks;
  public:
    PotentialMaster(Potential& p2, Box& box);
    virtual ~PotentialMaster() {}
    Box& getBox();
    virtual void computeAll(vector<PotentialCallback*> &callbacks);
    virtual void computeOne(int iAtom, double *ri, double &energy, bool isTrial);
    virtual void updateAtom(int iAtom) {}
    virtual void newAtom();
    virtual void removeAtom(int iAtom);
    double oldEnergy(int iAtom);
    void resetAtomDU();
    void processAtomU(int coeff);
    void addCallback(PotentialCallback* pcb);
};

class PotentialMasterCell : public PotentialMaster {
  protected:
    double range;
    int cellRange;
    double boxHalf[3];
    int numCells[3];
    vector<int> cellNextAtom;
    vector<int> atomCell;
    vector<int> cellLastAtom;
    int jump[3];
    vector<int> cellOffsets;
    vector<int> wrapMap;
    double** boxOffsets;
    int numAtoms;

    void handleComputeAll(int iAtom, int jAtom, double *ri, double *rj, double &ui, double &uj, double *fi, double *fj, double& uTot, double& virialTot, double rc2, bool doForces);
    void handleComputeOne(double *ri, double *rj, int jAtom, double& uTot, double rc2);
    int wrappedIndex(int i, int nc);
    void moveAtomIndex(int oldIndex, int newIndex);
  public:
    PotentialMasterCell(Potential& p2, Box& box, double potentialRange, int cellRange);
    ~PotentialMasterCell();
    void init();
    virtual void computeAll(vector<PotentialCallback*> &callbacks);
    virtual void computeOne(int iAtom, double *ri, double &energy, bool isTrial);
    virtual void updateAtom(int iAtom);
    virtual void newAtom();
    virtual void removeAtom(int iAtom);
    void assignCells();
    int cellForCoord(double *r);
    int* getNumCells();
};

class PotentialMasterList : public PotentialMasterCell {
  protected:
    int **nbrs;
    bool onlyUpNbrs; // standard MD only needs up.  MC or DMD needs down
    int *numAtomNbrsUp, *numAtomNbrsDn;
    int nbrsNumAtoms;
    int maxNab;
    bool forceReallocNbrs;
    double potentialRange;

    int checkNbrPair(int iAtom, int jAtom, double *ri, double *rj, double rc2);
  public:
    PotentialMasterList(Potential& p2, Box& box, double potentialRange, int cellRange, double nbrRange);
    ~PotentialMasterList();
    void reset();
    void setDoDownNbrs(bool doDown);

};
