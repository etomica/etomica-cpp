

#include "potential-master.h"
#include "bblock/system.h"
double PotentialMasterIntraMBnrg::computeOneMoleculeIntra(int iMolecule)
{
    bblock::System system;
    return system.OneBodyEnergy(false);
}


