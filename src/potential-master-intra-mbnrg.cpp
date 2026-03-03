

#include "potential-master.h"
#include "bblock/system.h"
#include "unit.h"

PotentialMasterIntraMBnrg::PotentialMasterIntraMBnrg(Box& b) : PotentialMasterIntra(b), pot("co2_archive")
{
        xyz1 = (double*)malloc(sizeof(double)*3*3);

}

PotentialMasterIntraMBnrg::~PotentialMasterIntraMBnrg()
{
        free(xyz1);
}

double PotentialMasterIntraMBnrg::computeOneMoleculeIntra(int iMolecule)
{
        int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
        box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
        for (int i = firstAtom; i <= lastAtom; i++)
        {
                double* p = box.getAtomPosition(i);
                xyz1[3*(i-firstAtom)] = p[0];
                xyz1[3*(i-firstAtom)+1] = p[1];
                xyz1[3*(i-firstAtom)+2] = p[2];
        }
        vector<double> energies = pot.eval(xyz1, 1);
        return Energy::toSim(energies[0]);

}




