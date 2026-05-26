#include "potential-molecular.h"
#include "bblock/system.h"
#include "unit.h"
#include "vector.h"

PotentialMolecularMBnrg3body::PotentialMolecularMBnrg3body() : pot("co2", "co2", "co2"),
                                                               grad(9, 0), chg_2b(6, 0), chggrad_2b(6, 0), pol_2b(6, 0), polfac_2b(6, 0),
                                                               chg(9, 0), chggrad(9, 0), pol(9, 0), polfac(9, 0)
{
    xyz1 = (double*)malloc(sizeof(double) * 3 * 3);
    xyz2 = (double*)malloc(sizeof(double) * 3 * 3);
    xyz3 = (double*)malloc(sizeof(double) * 3 * 3);
    {
        vector<double> sys_xyz(9*3, 0.0);
        nlohmann::json monomers_j_{};
        vector<size_t> first_index = {0, 0, 0};
        systools::SetCharges(sys_xyz, chg, "co2_archive", 3, 3, 0, chggrad, monomers_j_);
        systools::SetPolfac(polfac, "co2_archive", 3, 3, 0, monomers_j_);
        systools::SetPol(pol, "co2_archive", 3, 3, 0, monomers_j_);
        vector<size_t> sites = {3, 3, 3};
        vector<string> mon_id = {"co2_archive", "co2_archive", "co2_archive"};
        vector<size_t> num_atoms = {3, 3, 3};
        vector<pair<string, size_t>> mon_type_count = {{"co2_archive", 3}};
        vector<size_t> islocal = {1, 1, 1};
        vector<int> tags = {0, 1, 2};
        vector<double> edge_vectors;
        electrostatics.Initialize(chg, chggrad, polfac, pol, sys_xyz, mon_id, sites, first_index, mon_type_count,
                                   islocal, tags, false, 1E-16, 100, "iter", edge_vectors);

    }

    {
        vector<double> sys_xyz(6*3, 0.0);
        nlohmann::json monomers_j_{};
        vector<size_t> first_index = {0, 0};
        systools::SetCharges(sys_xyz, chg_2b, "co2_archive", 2, 3, 0, chggrad_2b, monomers_j_);
        systools::SetPolfac(polfac_2b, "co2_archive", 2, 3, 0, monomers_j_);
        systools::SetPol(pol_2b, "co2_archive", 2, 3, 0, monomers_j_);
        vector<size_t> sites = {3, 3};
        vector<string> mon_id = {"co2_archive", "co2_archive"};
        vector<size_t> num_atoms = {3, 3};
        vector<pair<string, size_t>> mon_type_count = {{"co2_archive", 2}};
        vector<size_t> islocal = {1, 1};
        vector<int> tags = {0, 1};
        vector<double> edge_vectors;


        electrostatics_pair.Initialize(chg_2b, chggrad_2b, polfac_2b, pol_2b, sys_xyz, mon_id, sites, first_index, mon_type_count,
                                   islocal, tags, false, 1E-16, 100, "iter", edge_vectors);

    }

}

PotentialMolecularMBnrg3body::~PotentialMolecularMBnrg3body()
{
    free(xyz1);
    free(xyz2);
    free(xyz3);

}

double PotentialMolecularMBnrg3body::u(Box& box, int iMolecule, int jMolecule, int kMolecule)
{
    // if (true) return 0;
    vector<double> xyz(9*3, 0.0);
    vector<double> xyz_2b(6*3, 0.0);

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
    int kSpecies, kMoleculeInSpecies;
    box.getMoleculeInfo(kMolecule, kSpecies, kMoleculeInSpecies, firstAtom, lastAtom);
    for (int k = firstAtom; k <= lastAtom; k++)
    {
        double* p = box.getAtomPosition(k);
        xyz3[3*(k-firstAtom)] = p[0];
        xyz3[3*(k-firstAtom)+1] = p[1];
        xyz3[3*(k-firstAtom)+2] = p[2];
        xyz[3*(k-firstAtom)+18] = p[0];
        xyz[3*(k-firstAtom)+19] = p[1];
        xyz[3*(k-firstAtom)+20] = p[2];

    }

    double energy3b = pot.eval(xyz1, xyz2, xyz3, 1);
    double triplet_total_energy = energy3b;
    vector<double> edge_vectors;
    electrostatics.SetNewParameters(xyz, chg, chggrad, pol, polfac, "iter", false, edge_vectors, std::numeric_limits<double>::infinity());
    triplet_total_energy += electrostatics.GetElectrostatics(grad, nullptr, false);
    // printf("triplet total energy: %f \n", triplet_total_energy);
    double dr23[3];
    for (int k=0; k<3; k++) dr23[k] = xyz3[k]-xyz2[+k];
    double r23 = Vector::squared(dr23);
    double dr12[3];
    double dr13[3];
    for (int k=0; k<3; k++) dr12[k] = xyz2[k]-xyz1[k];
    double r12 = Vector::squared(dr12);
    for (int k=0; k<3; k++) dr13[k] = xyz3[k]-xyz1[k];
    double r13 = Vector::squared(dr13);

    double minDist = sqrt(min({r12, r13, r23}));

    if (minDist > 40) return Energy::toSim(energy3b);
    // if (true) return Energy::toSim(energy3b);
    int molecules[3] = {iMolecule, jMolecule, kMolecule};
    for (int i = 0; i <= 2; i++)
    {
        int iSpecies, iMoleculeInSpecies, firstAtom, lastAtom;
        box.getMoleculeInfo(molecules[i], iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
        for (int z = firstAtom; z <= lastAtom; z++)
        {
            double* p = box.getAtomPosition(z);
            xyz_2b[3*(z-firstAtom)] = p[0];
            xyz_2b[3*(z-firstAtom)+1] = p[1];
            xyz_2b[3*(z-firstAtom)+2] = p[2];
        }

        for (int j = i + 1; j <= 2; j++)
        {
            box.getMoleculeInfo(molecules[j], iSpecies, iMoleculeInSpecies, firstAtom, lastAtom);
            for (int y = firstAtom; y <= lastAtom; y++)
            {
                double* p = box.getAtomPosition(y);
                xyz_2b[9+3*(y-firstAtom)] = p[0];
                xyz_2b[3*(y-firstAtom)+10] = p[1];
                xyz_2b[3*(y-firstAtom)+11] = p[2];

            }

            electrostatics_pair.SetNewParameters(xyz_2b, chg_2b, chggrad_2b, pol_2b, polfac_2b, "iter", false, edge_vectors, std::numeric_limits<double>::infinity());
            triplet_total_energy -= electrostatics.GetElectrostatics(grad, nullptr, false);


        }
    }


    return Energy::toSim(triplet_total_energy);
}
