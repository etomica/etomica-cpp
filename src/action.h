/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class Box;

class ConfigurationLattice {
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

class ConfigurationFile {
  protected:
    Box& box;
    const char* filename;

  public:
    ConfigurationFile(Box& box, const char* filename);
    ~ConfigurationFile();

    void go();
};

class Replicate {
  public:
    static void go(Box& box, const int replicates[3]);
};

class WriteXYZ {
  public:
    static void go(const char* filename, Box& box, char** symbols, bool append);
};
