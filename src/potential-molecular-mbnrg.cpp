#include "potential-molecular.h"
#include "bblock/system.h"
#include "unit.h"

PotentialMolecularMBnrg::PotentialMolecularMBnrg() : pot("co2", "co2"), grad(9, 0), chg(6, 0), chggrad(6, 0), pol(6, 0), polfac(6, 0)
{
    xyz1 = (double*)malloc(sizeof(double) * 3 * 3);
    xyz2 = (double*)malloc(sizeof(double) * 3 * 3);
    vector<double> sys_xyz(6*3, 0.0);
    nlohmann::json monomers_j_{};
    vector<size_t> first_index = {0, 0};
    systools::SetCharges(sys_xyz, chg, "co2", 2, 3, 0, chggrad, monomers_j_);
    systools::SetPolfac(polfac, "co2", 2, 3, 0, monomers_j_);
    systools::SetPol(pol, "co2", 2, 3, 0, monomers_j_);
    vector<size_t> sites = {3, 3};
    vector<string> mon_id = {"co2", "co2"};
    vector<size_t> num_atoms = {3, 3};
    vector<pair<string, size_t>> mon_type_count = {{"co2", 2}};
    vector<size_t> islocal = {0, 0};
    vector<int> tags = {0, 1};
    vector<double> edge_vectors;
    vector<double> c6_lr(6, 0);

    electrostatics.Initialize(chg, chggrad, polfac, pol, sys_xyz, mon_id, sites, first_index, mon_type_count,
                               islocal, tags, false, 1E-16, 100, "iter", edge_vectors);
    systools::SetC6LongRange(c6_lr, "co2", 2, 3, 0, monomers_j_);


    dispersion.Initialize(c6_lr, sys_xyz, mon_id, num_atoms, mon_type_count, islocal, false, edge_vectors);

}

PotentialMolecularMBnrg::~PotentialMolecularMBnrg()
{
    free(xyz1);
    free(xyz2);
}

double PotentialMolecularMBnrg::u(Box& box, int iMolecule, int jMolecule)
{
    vector<double> xyz(6*3, 0.0);
    int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
    for (int i = firstAtom; i <= lastAtom; i++)
    {
        double* p = box.getAtomPosition(i);
        xyz1[3*(i-firstAtom)] = p[0];
        xyz1[3*(i-firstAtom)+1] = p[1];
        xyz1[3*(i-firstAtom)+2] = p[2];
        xyz[3*(i-firstAtom)] = p[0];
        xyz[3*(i-firstAtom)+1] = p[1];
        xyz[3*(i-firstAtom)+2] = p[2];
    }

    int jSpecies, jMoleculeInSpecies;
    box.getMoleculeInfo(jMolecule, jSpecies, jMoleculeInSpecies, firstAtom, lastAtom);
    for (int j = firstAtom; j <= lastAtom; j++)
    {
        double* p = box.getAtomPosition(j);
        xyz2[3*(j-firstAtom)] = p[0];
        xyz2[3*(j-firstAtom)+1] = p[1];
        xyz2[3*(j-firstAtom)+2] = p[2];
        xyz[9+3*(j-firstAtom)] = p[0];
        xyz[3*(j-firstAtom)+10] = p[1];
        xyz[3*(j-firstAtom)+11] = p[2];

    }
    double energy2b = pot.eval(xyz1, xyz2, 1);
    double total_energy = energy2b;
    if (energy2b < 6000)
    {
        vector<double> edge_vectors;
        vector<pair<string, string>> ignore_disp;
        dispersion.SetNewParameters(xyz, ignore_disp, false, std::numeric_limits<double>::infinity(), edge_vectors);
        double disp_energy = dispersion.GetDispersion(grad, nullptr, false);
        electrostatics.SetNewParameters(xyz, chg, chggrad, pol, polfac, "iter", false, edge_vectors, std::numeric_limits<double>::infinity());
        double elec_energy = electrostatics.GetElectrostatics(grad, nullptr, false);
        total_energy += elec_energy + disp_energy;
        printf("%f %f %f \n", energy2b, elec_energy, disp_energy);

    }
    else
    {
        printf("%f \n", energy2b);

    }
    return Energy::toSim(total_energy);
}
