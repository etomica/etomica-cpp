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
    PotentialCallbackPressure(Box& box, double temperature, bool takesForces);
    ~PotentialCallbackPressure() {}
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};
    
class PotentialCallbackEnergy : public PotentialCallback {
  protected:
    double data[1];
  public:
    PotentialCallbackEnergy();
    ~PotentialCallbackEnergy() {}
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
    ~PotentialCallbackHMA();
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
    double temperature;
    double Pharm;
    double* data;
    double** latticePositions;
    double** latticeOrientations;
    bool returnAnh, computingLat;
    double uLat, pLat;
  public:
    PotentialCallbackMoleculeHMA(Box& box, SpeciesList& speciesList, double temperature, double Pharm);
    ~PotentialCallbackMoleculeHMA();
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
    void setReturnAnharmonic(bool returnAnharmonic, PotentialMaster* potentialMaster);
};

