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
    double energy2b = pot.eval(xyz1, xyz2, 1);
    monomers_j_ = {};
    vector<double> sys_xyz;
    sys_xyz.insert(sys_xyz.end(), xyz1[0], xyz1[8]);
    sys_xyz.insert(sys_xyz.end(), xyz2[0], xyz2[8]);
    systools::SetC6LongRange(c6_lr, "co2", 2, 3, 0, monomers_j_);
    disp::Dispersion dispersion;

    vector<string> mon_id = {"co2", "co2"};
    vector<size_t> num_atoms = {3, 3};
    vector<pair<string, size_t>> mon_type_count = {{"co2", 2}};
    vector<size_t> islocal = {0, 0};
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            edge_vectors.push_back(box.getEdgeVector(i)[j]);
        }
    }

    dispersion.Initialize(c6_lr, sys_xyz, mon_id, num_atoms, mon_type_count, islocal, false, edge_vectors);
    double disp_energy = dispersion.GetDispersion(grad, virial, false);


    // lj::LennardJones lj;
    // systools::SetLJLongRange(lj_lr, "co2", 2, 3, 0, );
    // lj.Initialize(lj_lr, sys_xyz, mon_id, num_atoms, mon_type_count, islocal, false, edge_vectors);
    // double lj_energy = lj.GetLennardJones(grad, virial, false);


    elec::Electrostatics electrostatics;
    vector<size_t> first_index = {0, 0};
    systools::SetCharges(sys_xyz, chg, "co2", 2, 3, 0, chggrad, monomers_j_);
    systools::SetPolfac(polfac, "co2", 2, 3, 0, monomers_j_);
    systools::SetPol(pol, "co2", 2, 3, 0, monomers_j_);
    vector<size_t> sites = {3, 3};
    electrostatics.Initialize(chg, chggrad, polfac, pol, sys_xyz, mon_id, sites, first_index, mon_type_count,
                               islocal, tags, false, 1E-16, 100, "iter", edge_vectors);
    double elec_energy = electrostatics.GetElectrostatics(grad, virial, false);
    int kSpecies, kMoleculeInSpecies;
    box.getMoleculeInfo(kMolecule, kSpecies, kMoleculeInSpecies, firstAtom, lastAtom);
    for (int k = firstAtom; k <= lastAtom; k++)
    {
        double* p = box.getAtomPosition(k);
        xyz3[3*(k-firstAtom)] = p[0];
        xyz3[3*(k-firstAtom)+1] = p[1];
        xyz3[3*(k-firstAtom)+2] = p[2];
    }
    double energy3b = pot3.eval(xyz1, xyz2, xyz3, 3);


    double total_energy = energy2b + elec_energy + disp_energy + energy3b;
    return Energy::toSim(total_energy);
}
