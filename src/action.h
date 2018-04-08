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

class ConfigurationFile : Action {
  protected:
    Box& box;
    const char* filename;

  public:
    ConfigurationFile(Box& box, const char* filename);
    ~ConfigurationFile();
    virtual void go();
};
