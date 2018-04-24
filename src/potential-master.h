#pragma once

#include <vector>
#include <math.h>
#include <cstddef>
#include <set>
#include <algorithm>
#include <complex>
#include "box.h"
#include "potential.h"
#include "potential-angle.h"
#include "potential-molecular.h"
#include "potential-callback.h"

using namespace std;

class EmbedF {
  public:
    EmbedF() {}
    virtual ~EmbedF() {}
    virtual double f(double rhoSum) = 0;
    virtual void f012(double rhoSum, double &f, double &df, double &d2f) = 0;
};

class EmbedFsqrt : public EmbedF {
  private:
    const double eps;
  public:
    EmbedFsqrt(double Ceps) : EmbedF(), eps(Ceps) {}
    ~EmbedFsqrt() {}
    double f(double rhoSum) { return -eps*sqrt(rhoSum); }
    void f012(double rhoSum, double &f, double &df, double &d2f) {
      double s = sqrt(rhoSum);
      f = -eps*s;
      df = -0.5*eps/s;
      d2f = 0;
    }
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
    int numRawBoxOffsets;
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
    Potential** rhoPotentials;
    EmbedF **embedF;
    int* numAtomsByType;
    double** pairCutoffs;
    double *rhoCutoffs;
    Box& box;
    vector<double> uAtom;
    vector<double> duAtom;
    bool duAtomSingle, duAtomMulti;
    vector<int> uAtomsChanged;
    double** force;
    int numForceAtoms, numRhoSumAtoms;
    double* rhoSum;
    double* idf;
    vector<double> rdrho;
    vector<int> rhoAtomsChanged;
    vector<double> drhoSum;

    vector<PotentialCallback*> pairCallbacks;
    const int numAtomTypes;
    // one vector<vector<int*>> for each species
    // each species has a list of bonded pairs for each potential
    vector<vector<int*> > *bondedPairs;
    vector<int> **bondedAtoms;
    vector<Potential*> *bondedPotentials;
    vector<vector<int*> > *bondAngleTriplets;
    vector<PotentialAngle*> *bondAnglePotentials;
    const bool pureAtoms;
    bool rigidMolecules;
    bool doTruncationCorrection, doSingleTruncationCorrection;
    const bool embeddingPotentials;
    double* charges;
    double kBasis[3];
    double kCut, alpha;
    complex<double> *sFacAtom;
    vector<complex<double>> sFac;
    vector<complex<double>> eik[3];
    vector<complex<double>> dsFacMolecule;
    vector<double> fExp;
    bool doEwald;
    double minR2;

    void computeOneMoleculeBonds(const int iSpecies, const int iMolecule, double &u1);
    void handleOneBondPair(bool doForces, double &uTot, int iAtom, int jAtom, Potential* p);
    void handleOneBondAngleTriplet(bool doForces, double &uTot, int iAtom, int jAtom, int kATom, PotentialAngle* p);
    void computeAllBonds(bool doForces, double &uTot);
    void computeAllTruncationCorrection(double &uTot, double &virialTot);
    double computeOneTruncationCorrection(const int iAtom);
    inline bool checkSkip(int jAtom, int iSpecies, int iMolecule, vector<int> *iBondedAtoms) {
      if (pureAtoms) return false;
      int jMolecule, jFirstAtom, jSpecies;
      box.getMoleculeInfoAtom(jAtom, jMolecule, jSpecies, jFirstAtom);
      if (rigidMolecules) return iSpecies==jSpecies && iMolecule == jMolecule;
      return binary_search(iBondedAtoms->begin(), iBondedAtoms->end(), jAtom-jFirstAtom);
    }
    void handleComputeAll(int iAtom, int jAtom, const double *ri, const double *rj, const double *jbo, Potential* pij, double &ui, double &uj, double* fi, double* fj, double& uTot, double& virialTot, const double rc2, Potential* iRhoPotential, const double iRhoCutoff, const int iType, const int jType, const bool doForces, const bool skipIntra) {
      double dr[3];
      dr[0] = (rj[0]+jbo[0])-ri[0];
      dr[1] = (rj[1]+jbo[1])-ri[1];
      dr[2] = (rj[2]+jbo[2])-ri[2];
      double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      if (r2 < rc2 && (!skipIntra || r2 > minR2)) {
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
        if (doForces) {
          du /= r2;
          double x = dr[0]*du;
          fi[0] += x;
          fj[0] -= x;
          x = dr[1]*du;
          fi[1] += x;
          fj[1] -= x;
          x = dr[2]*du;
          fi[2] += x;
          fj[2] -= x;
        }
      }
      if (embeddingPotentials) {
        if (r2 < iRhoCutoff) {
          double rho, drho, d2rho;
          iRhoPotential->u012(r2, rho, drho, d2rho);
          rhoSum[iAtom] += rho;
          rdrho.push_back(drho);
          //if (iAtom==0||jAtom==0) printf("%d %d %f %f\n", iAtom, jAtom, sqrt(r2), drho);
          if (jType == iType) {
            rhoSum[jAtom] += rho;
          }
          else if (r2 < rhoCutoffs[jType]) {
            rhoPotentials[jType]->u012(r2, rho, drho, d2rho);
            rhoSum[jAtom] += rho;
            rdrho.push_back(drho);
          }
        }
        else if (r2 < rhoCutoffs[jType]) {
          double rho, drho, d2rho;
          rhoPotentials[jType]->u012(r2, rho, drho, d2rho);
          rhoSum[jAtom] += rho;
          rdrho.push_back(drho);
          //if (iAtom==0||jAtom==0) printf("%d %d %f %f\n", iAtom, jAtom, sqrt(r2), drho);
        }
      }
    }
    void handleComputeAllEmbed(const int iAtom, const int jAtom, const int iType, const int jType, const double *ri, const double *rj, const double *jbo, const double df, double &virialTot, Potential* iRhoPotential, double iRhoCutoff, int &rdrhoIdx) {
      double dr[3];
      dr[0] = (rj[0]+jbo[0])-ri[0];
      dr[1] = (rj[1]+jbo[1])-ri[1];
      dr[2] = (rj[2]+jbo[2])-ri[2];
      double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

      if (r2 < iRhoCutoff) {
        // if rij < cutoff for rhoi, then we have a force
        if (iType==jType) {
          double fac = (df + idf[jAtom]) * rdrho[rdrhoIdx];
          rdrhoIdx++;
          virialTot += fac;
          fac /= r2;
          //if (iAtom==0||jAtom==0) printf("go %d %d  %f %f  %f %f  %f\n", iAtom, jAtom, sqrt(r2), rdrho[rdrhoIdx], df, idf[jAtom], fac);
          for (int k=0; k<3; k++) {
            double fk = dr[k] * fac;
            force[iAtom][k] += fk;
            force[jAtom][k] -= fk;
          }
        }
        else {
          double fac = df * rdrho[rdrhoIdx];
          rdrhoIdx++;
          if (r2 < rhoCutoffs[jAtom]) {
            fac += idf[jAtom] * rdrho[rdrhoIdx];
            rdrhoIdx++;
          }
          virialTot += fac;
          fac /= r2;
          for (int k=0; k<3; k++) {
            double fk = dr[k] * fac;
            force[iAtom][k] += fk;
            force[jAtom][k] -= fk;
          }
        }
      }
      else if (r2 < rhoCutoffs[jType]) {
        // if rij < cutoff for rhoi, then we have a force
        double fac = idf[jAtom] * rdrho[rdrhoIdx];
        rdrhoIdx++;
        virialTot += fac;
        fac /= r2;
        for (int k=0; k<3; k++) {
          double fk = dr[k] * fac;
          force[iAtom][k] += fk;
          force[jAtom][k] -= fk;
        }
      }
    }
    void handleOldEmbedding(const double *ri, const double *rj, const double *jbo, const int jAtom, double& uTot, const int jType) {
      double dx = ri[0]-(rj[0]+jbo[0]);
      double dy = ri[1]-(rj[1]+jbo[1]);
      double dz = ri[2]-(rj[2]+jbo[2]);
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < rhoCutoffs[jType]) {
        rhoAtomsChanged.push_back(jAtom);
        double rho = rhoPotentials[jType]->u(r2);
        drhoSum[jAtom] = -rho;
        // we need to compute the old energy of the system with and without iAtom
        uTot -= embedF[jType]->f(rhoSum[jAtom]-rho);
        uTot += embedF[jType]->f(rhoSum[jAtom]);
        //printf("%d %d %f  %e %e   %e %e  %e\n", iAtom, jAtom, r2, rhoSum[jAtom], rho, embedF[jType]->f(rhoSum[jAtom]), embedF[jType]->f(rhoSum[jAtom]-rho), embedF[jType]->f(rhoSum[jAtom])-embedF[jType]->f(rhoSum[jAtom]-rho));
      }
    }
    void handleComputeOne(Potential* pij, const double *ri, const double *rj, const double* jbo, const int iAtom, const int jAtom, double& uTot, double rc2, const double iRhoCutoff, Potential* iRhoPotential, const int iType, const int jType, const bool skipIntra) {
      double dx = ri[0]-(rj[0]+jbo[0]);
      double dy = ri[1]-(rj[1]+jbo[1]);
      double dz = ri[2]-(rj[2]+jbo[2]);
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < rc2 && (!skipIntra || r2 > minR2)) {
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
      if (embeddingPotentials) {
        if (r2 < iRhoCutoff) {
          double rho = iRhoPotential->u(r2);
          drhoSum[iAtom] += rho;
          if (iType == jType) {
            if (drhoSum[jAtom] == 0) {
              rhoAtomsChanged.push_back(jAtom);
            }
            // we need the energy of the configuration with the atom in the new spot
            // minus the energy with no atom
            uTot += embedF[jType]->f(rhoSum[jAtom]+drhoSum[jAtom]+rho);
            uTot -= embedF[jType]->f(rhoSum[jAtom]+drhoSum[jAtom]);
            // drhoSum will hold new-old
            // processAtom(+1) handles this and we should be done
            // we'll be called again and recompute rho for the old config
            // for processAtom(-1), we ignore drhoSum since it was already handled
            drhoSum[jAtom] += rho;
          }
          else if (r2 < rhoCutoffs[jType]) {
            double rho = rhoPotentials[jType]->u(r2);
            if (drhoSum[jAtom] == 0) rhoAtomsChanged.push_back(jAtom);
            uTot += embedF[jType]->f(rhoSum[jAtom]+drhoSum[jAtom]+rho);
            uTot -= embedF[jType]->f(rhoSum[jAtom]+drhoSum[jAtom]);
            drhoSum[jAtom] += rho;
          }
        }
        else if (r2 < rhoCutoffs[jType]) {
          double rho = rhoPotentials[jType]->u(r2);
          if (drhoSum[jAtom] == 0) rhoAtomsChanged.push_back(jAtom);
          uTot += embedF[jType]->f(rhoSum[jAtom]+drhoSum[jAtom]+rho);
          uTot -= embedF[jType]->f(rhoSum[jAtom]+drhoSum[jAtom]);
          drhoSum[jAtom] += rho;
        }
      }
    }
    virtual void computeOneInternal(const int iAtom, const double *ri, double &u1, const int iSpecies, const int iMolecule, const int iFirstAtom, const bool onlyAtom);
    virtual double oldEmbeddingEnergy(int iAtom);
    void computeAllFourier(const bool doForces, double &uTot);
    double oneMoleculeFourierEnergy(int iMolecule, bool oldEnergy);
    double computeFourierIntramolecular(int iMolecule, bool doForces);

  public:
    PotentialMaster(const SpeciesList &speciesList, Box& box, bool doEmbed);
    virtual ~PotentialMaster();
    Box& getBox();
    void setDoTruncationCorrection(bool doCorrection);
    void setDoSingleTruncationCorrection(bool doCorrection);
    virtual void setPairPotential(int iType, int jType, Potential* pij);
    virtual void setRhoPotential(int jType, Potential* rhoj);
    virtual void setEmbedF(int iType, EmbedF* Fi);
    void setBondPotential(int iSpecies, vector<int*> &bondedPairs, Potential *pBond);
    void setBondAnglePotential(int iSpecies, vector<int*> &bondedTriplets, PotentialAngle *pBondAngle);
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
    virtual double oldIntraMoleculeEnergyLS(int iAtom, int iLastAtom) {return 0;}
    void resetAtomDU();
    void processAtomU(int coeff);
    void addCallback(PotentialCallback* pcb);
    virtual double uTotalFromAtoms();
    void setCharge(int iType, double charge);
    void setEwald(double kCut, double alpha);
    virtual void updateVolume() {}
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

    virtual void computeOneInternal(const int iAtom, const double *ri, double &energy, const int iSpecies, const int iMolecule, const int iFirstAtom, const bool onlyAtom);
    virtual double oldEmbeddingEnergy(int iAtom);

  public:
    PotentialMasterCell(const SpeciesList &speciesList, Box& box, bool doEmbed, int cellRange);
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
    virtual double oldIntraMoleculeEnergyLS(int iAtom, int iLastAtom);
    int* getNumCells();
    virtual void updateVolume();
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

    int checkNbrPair(int iAtom, int jAtom, const bool skipIntra, double *ri, double *rj, double rc2, double minR2, double *jbo);
  public:
    PotentialMasterList(const SpeciesList& speciesList, Box& box, bool doEmbed, int cellRange, double nbrRange);
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

