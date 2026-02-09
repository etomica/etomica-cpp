/* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <algorithm>
#include <stdio.h>
#include "potential-master.h"

#include <limits>

#include "potential-angle.h"
#include "alloc2d.h"
#include "util.h"

PotentialMasterIntra::PotentialMasterIntra(Box& box) : box(box)
{
}

void PotentialMasterIntra::setMoleculeEnergy(int iMolecule, double u1)
{
    uMolecule[iMolecule] = u1;
}
double PotentialMasterIntra::oldMoleculeEnergyIntra(int iMolecule)
{
    // We may not have saved energies at the beginning of a simulation, so need to check for nan
    if (!std::isnan(uMolecule[iMolecule]))
    {
        // double uCheck = computeOneMoleculeIntra(iMolecule);
        // if (fabs(uCheck - uMolecule[iMolecule])>0.001)
        // {
        //   printf("%f %f \n", uCheck, uMolecule[iMolecule]);
        //   abort();
        // }
        return uMolecule[iMolecule];
    }
    double u1 = computeOneMoleculeIntra(iMolecule);
    uMolecule[iMolecule] = u1;
    return u1;

}
void PotentialMasterIntra::init() {
    int nm = box.getTotalNumMolecules();
    int s = uMolecule.size();
    if (s < nm || s > 2*nm)
    {
        uMolecule.resize(nm);
    }
    std::fill(uMolecule.begin(), uMolecule.end(), NAN);
}
