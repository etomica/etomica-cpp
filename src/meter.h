#pragma once

#include "average.h"
#include "integrator.h"

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
    IntegratorMC& integrator;
    double data[1];
  public:
    MeterPotentialEnergy(IntegratorMC& integrator);
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

class MeterFullCompute : public Meter {
  protected:
    PotentialMaster& potentialMaster;
    double *data;
    vector<PotentialCallback*> callbacks;
  public:
    MeterFullCompute(PotentialMaster& potentialMaster);
    ~MeterFullCompute() {}
    double* getData();
    void addCallback(PotentialCallback* pcb);
};

class DataPump : public IntegratorListener {
  private:
    Meter& meter;
    int interval;
    int intervalCountdown;
    DataSink* sink;
  public:
    DataPump(Meter& meter, int interval);
    DataPump(Meter& meter, int interval, DataSink* sink);
    ~DataPump();
    void stepFinished();
    DataSink* getDataSink();
};
