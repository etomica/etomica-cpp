#include "potential-molecular.h"
#include "bblock/system.h"
#include "unit.h"

PotentialMolecularMBnrg3body::PotentialMolecularMBnrg3body() : pot("co2", "co2", "co2")
{
    xyz1 = (double*)malloc(sizeof(double) * 3 * 3);
    xyz2 = (double*)malloc(sizeof(double) * 3 * 3);
    xyz3 = (double*)malloc(sizeof(double) * 3 * 3);

}

PotentialMolecularMBnrg3body::~PotentialMolecularMBnrg3body()
{
    free(xyz1);
    free(xyz2);
    free(xyz3);

}

double PotentialMolecularMBnrg3body::u(Box& box, int iMolecule, int jMolecule, int kMolecule)
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
    int jSpecies, jMoleculeInSpecies;
    box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, firstAtom, lastAtom);
    for (int j = firstAtom; j <= lastAtom; j++)
    {
        double* p = box.getAtomPosition(j);
        xyz2[3*(j-firstAtom)] = p[0];
        xyz2[3*(j-firstAtom)+1] = p[1];
        xyz2[3*(j-firstAtom)+2] = p[2];
    }
    int kSpecies, kMoleculeInSpecies;
    box.getMoleculeInfo(kMolecule, kSpecies, kMoleculeInSpecies, firstAtom, lastAtom);
    for (int k = firstAtom; k <= lastAtom; k++)
    {
        double* p = box.getAtomPosition(k);
        xyz3[3*(k-firstAtom)] = p[0];
        xyz3[3*(k-firstAtom)+1] = p[1];
        xyz3[3*(k-firstAtom)+2] = p[2];
    }
    double energy3b = pot.eval(xyz1, xyz2, xyz3, 3);
    return Energy::toSim(energy3b);
}
