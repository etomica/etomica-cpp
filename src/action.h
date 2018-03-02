#pragma once

class Box;

class Action {
  public:
    Action() {}
    virtual ~Action() {}
    virtual void go() = 0;
};

class ConfigurationLattice : Action {
  protected:
    Box& box;
    int numBasisAtoms;
    double** basis;
    double* cellShape;

  public:
    ConfigurationLattice(Box& box, double** basis, double* cellShape);
    ~ConfigurationLattice();
    virtual void go();
};
