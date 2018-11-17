/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class Box;
class SpeciesList;
class PotentialMaster;

class PotentialCallback {
  public:
    bool callPair;
    bool callFinished;
    bool takesForces;
    bool takesPhi;
    bool takesDFDV;

    PotentialCallback();
    virtual ~PotentialCallback() {}
    virtual void reset() {}
    virtual void pairCompute(int iAtom, int jAtom, double* dr, double u, double du, double d2u) {}
    virtual void pairComputePhi(int iAtom, int jAtom, double phi[3][3]) {}
    virtual void computeDFDV(int iAtom, double* dFdV) {}
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
    PotentialCallbackPressure(Box& box, double temperature, bool takesForces);
    virtual ~PotentialCallbackPressure() {}
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};

class PotentialCallbackEnergy : public PotentialCallback {
  protected:
    double data[1];
  public:
    PotentialCallbackEnergy();
    virtual ~PotentialCallbackEnergy() {}
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};

class PotentialCallbackHMA : public PotentialCallback {
  protected:
    Box& box;
    double temperature;
    double Pharm;
    double* data;
    double** latticePositions;
    double phiSum;
    bool doD2;
    double dr0[3];
    bool returnAnh, computingLat;
    double uLat, pLat;
  public:
    PotentialCallbackHMA(Box& box, double temperature, double Pharm, bool doD2);
    virtual ~PotentialCallbackHMA();
    virtual void pairCompute(int iAtom, int jAtom, double* dr, double u, double du, double d2u);
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
    void setReturnAnharmonic(bool returnAnharmonic, PotentialMaster* potentialMaster);
};

class PotentialCallbackMoleculeHMA : public PotentialCallback {
  protected:
    Box& box;
    SpeciesList& speciesList;
    PotentialMaster* potentialMaster;
    double temperature;
    double Pharm;
    double* data;
    double** latticePositions;
    double** latticeOrientations;
    double** dFdV;
    double** phiTotal;
    double** com;
    double** dRdV;
    bool returnAnh, computingLat, computingPshift;
    double uLat, pLat;

    void computeShift(double** f);

  public:
    PotentialCallbackMoleculeHMA(Box& box, SpeciesList& speciesList, PotentialMaster* pm, double temperature, double Pharm);
    virtual ~PotentialCallbackMoleculeHMA();
    virtual void findShiftV();
    virtual double** getDRDV();
    virtual void pairCompute(int iAtom, int jAtom, double* dr, double u, double du, double d2u);
    virtual void pairComputePhi(int iAtom, int jAtom, double phi[3][3]);
    virtual void computeDFDV(int iAtom, double* idFdV);
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
    void setReturnAnharmonic(bool returnAnharmonic);
};

class PotentialCallbackBulkModulus : public PotentialCallback {
  protected:
    Box& box;
    double temperature;
    double data[2];
  public:
    PotentialCallbackBulkModulus(Box& box, double temperature);
    virtual ~PotentialCallbackBulkModulus() {}
    virtual int getNumData();
    virtual void reset();
    virtual void pairCompute(int iAtom, int jAtom, double* dr, double u, double du, double d2u);
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual double* getData();
};

class PotentialCallbackMoleculePhi : public PotentialCallback {
  protected:
    Box& box;
    SpeciesList& speciesList;
    PotentialMaster& potentialMaster;
    int ndof;
    double** atomPhiTotal;
    double** moleculePhiTotal;
    int* nmap;
    double** com;
    double uLat, pLat;

  public:
    PotentialCallbackMoleculePhi(Box& box, SpeciesList& speciesList, PotentialMaster& pm);
    virtual ~PotentialCallbackMoleculePhi();
    virtual void reset();
    virtual void pairCompute(int iAtom, int jAtom, double* dr, double u, double du, double d2u);
    virtual void pairComputePhi(int iAtom, int jAtom, double phi[3][3]);
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    double** getMoleculePhi();
};
