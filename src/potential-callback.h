#pragma once

class Box;
class SpeciesList;

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
    double data[2];
    double** latticePositions;
  public:
    PotentialCallbackHMA(Box& box, double temperature, double Pharm);
    ~PotentialCallbackHMA();
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};

class PotentialCallbackMoleculeHMA : public PotentialCallback {
  protected:
    Box& box;
    SpeciesList& speciesList;
    double temperature;
    double data[1];
    double** latticePositions;
    double** latticeOrientations;
  public:
    PotentialCallbackMoleculeHMA(Box& box, SpeciesList& speciesList, double temperature);
    ~PotentialCallbackMoleculeHMA();
    virtual void allComputeFinished(double uTot, double virialTot, double** f);
    virtual int getNumData();
    virtual double* getData();
};

