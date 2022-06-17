/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "integrator.h"
#include "species.h"

class Integrator;
class IntegratorMD;

class Meter {
  protected:
    int nData;
  public:
    Meter(int n) : nData(n) {}
    virtual ~Meter() {}
    virtual double* getData() = 0;
    int getNumData() {return nData;}
};

class MeterPotentialEnergy : public Meter {
  private:
    Integrator& integrator;
    double data[1];
  public:
    MeterPotentialEnergy(Integrator& integrator);
    ~MeterPotentialEnergy() {}
    double* getData();
};

class MeterNumAtoms : public Meter {
  private:
    Box& box;
    double data[1];
  public:
    MeterNumAtoms(Box& box);
    ~MeterNumAtoms() {}
    double* getData();
};

class MeterDensity : public Meter {
  private:
    Box& box;
    double data[1];
  public:
    MeterDensity(Box& box);
    ~MeterDensity() {}
    double* getData();
};

class MeterKineticEnergy : public Meter {
  private:
    IntegratorMD* integrator;
    double data[2];
  public:
    MeterKineticEnergy(IntegratorMD* integrator);
    ~MeterKineticEnergy() {}
    double* getData();
};

class MeterFullCompute : public Meter {
  protected:
    PotentialMaster& potentialMaster;
    double *data;
    vector<PotentialCallback*> callbacks;
    bool doCompute;
  public:
    MeterFullCompute(PotentialMaster& potentialMaster);
    ~MeterFullCompute();
    double* getData();
    void addCallback(PotentialCallback* pcb);
    void setDoCompute(bool doCompute);
};

class MeterPressureFD : public Meter {
  protected:
    SpeciesList& speciesList;
    PotentialMaster& potentialMaster;
    double data[1];
    vector<PotentialCallback*> callbacks;
    double eps;
    double temperature;

    void scale(double s, double* boxSize);
  public:
    MeterPressureFD(SpeciesList& sl, PotentialMaster& potentialMaster, double temperature);
    ~MeterPressureFD();
    void setEps(double newEps);
    double* getData();
};

class MeterWidomInsertion : public Meter {
  protected:
    Box& box;
    int iSpecies;
    PotentialMaster& potentialMaster;
    Random& random;
    double data[1];
    double temperature;
    int numTrials;

  public:
    MeterWidomInsertion(Box& box, int iSpecies, PotentialMaster& potentialMaster, Random& rand, double temperature, int numTrials);
    ~MeterWidomInsertion();
    double* getData();
};

class MeterSteps : public Meter {
  private:
    Integrator& integrator;
    double data[1];
  public:
    MeterSteps(Integrator& integrator);
    ~MeterSteps() {}
    double* getData();
};

class MeterTime : public Meter {
  private:
    IntegratorMD& integrator;
    double data[1];
  public:
    MeterTime(IntegratorMD& integrator);
    ~MeterTime() {}
    double* getData();
};

class MeterRDF : public Meter {
  private:
    Box& box;
    double rMax;
    double *data;
    double *rData;
  public:
    MeterRDF(Box& box, double rMax, int nBins);
    virtual ~MeterRDF();
    void dispose();
    void reset();
    void setup();
    double* getData();
    double* getRData();
    void setNBins(int nBins);
    void setRMax(double rMax);
};
