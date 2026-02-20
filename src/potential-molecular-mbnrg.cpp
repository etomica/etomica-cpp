#include "potential-molecular.h"
#include "bblock/system.h"
#include "unit.h"

PotentialMolecularMBnrg::PotentialMolecularMBnrg() : pot("co2", "co2")
{
    xyz1 = (double*)malloc(sizeof(double) * 3 * 3);
    xyz2 = (double*)malloc(sizeof(double) * 3 * 3);
}

PotentialMolecularMBnrg::~PotentialMolecularMBnrg()
{
    free(xyz1);
    free(xyz2);
}

double PotentialMolecularMBnrg::u(Box& box, int iMolecule, int jMolecule)
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
    double energy = pot.eval(xyz1, xyz2, 1);
    printf("%f %f %f\n", xyz1[0], xyz1[1], xyz1[2]);
    printf("%f %f %f\n", xyz1[3], xyz1[4], xyz1[5]);
    printf("%f %f %f\n", xyz1[6], xyz1[7], xyz1[8]);
    printf("%f %f %f\n", xyz2[0], xyz2[1], xyz2[2]);
    printf("%f %f %f\n", xyz2[3], xyz2[4], xyz2[5]);
    printf("%f %f %f\n", xyz2[6], xyz2[7], xyz2[8]);

    printf("%f\n", Energy::toSim(energy));

    return Energy::toSim(energy);
}
