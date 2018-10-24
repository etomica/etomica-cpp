/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class AtomInfo {
  private:
    int numAtomTypes;
    double *mass;
    
  public:
    AtomInfo();
    virtual ~AtomInfo();
    int addAtomType(double mass);
    double getMass(int iType) const;
    int getNumTypes() const;
};

