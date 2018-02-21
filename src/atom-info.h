#pragma once

class AtomInfo {
  private:
    int numAtomTypes;
    double *mass;
    
  public:
    AtomInfo();
    virtual ~AtomInfo();
    int addAtomType(double mass);
    double getMass(int iType);
};

