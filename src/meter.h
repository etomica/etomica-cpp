/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "data-sink.h"
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
    MeterKineticEnergy();
    ~MeterKineticEnergy() {}
    void setIntegrator(IntegratorMD* integrator);
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
    SpeciesList speciesList;
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
