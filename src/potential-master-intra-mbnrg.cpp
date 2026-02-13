

#include "potential-master.h"
#include "bblock/system.h"

double PotentialMasterIntraMBnrg::computeOneMoleculeIntra(int iMolecule)
{
    const int n_monomers = 1;                                                                                 \
    std::vector<double> real_coords{-3.7350600000e-03, 2.1493800000e-02,  -2.0987560000e-02,                  \
                                    1.1921700000e-03,  -7.0306800000e-03, 1.2319206500e+00,                   \
                                    1.1921700000e-03,  -7.0306800000e-03, -1.2182866600e+00};                 \
    std::vector<std::string> atom_names{"C", "O", "O"};                                                       \
    std::vector<std::string> monomer_names = {"co2_archive"};                                                 \
    std::vector<size_t> n_atoms_vector = {3};                                                                 \
    eff::Bond co2_bond1("bond", std::vector<size_t>{1, 2}, "morse");                                          \
    eff::Bond co2_bond2("bond", std::vector<size_t>{1, 3}, "morse");                                          \
    eff::Angles co2_angle1("angle", std::vector<size_t>{2, 1, 3}, "harm");                                    \
    std::vector<double> bond_morse_linear_parameters = {189.05};                                              \
    std::vector<double> bond_morse_nonlinear_parameters = {2.57, 1.16};                                       \
    std::vector<double> angle_harm_linear_parameters1 = {87.87};                                              \
    std::vector<double> angle_harm_nonlinear_parameters1 = {3.14159265};
    co2_bond1.SetParameters(bond_morse_linear_parameters, bond_morse_nonlinear_parameters);
    co2_bond2.SetParameters(bond_morse_linear_parameters, bond_morse_nonlinear_parameters);
    std::vector<eff::Bond> bond_vec = {co2_bond1, co2_bond2};
    co2_angle1.SetParameters(angle_harm_linear_parameters1, angle_harm_nonlinear_parameters1);
    std::vector<eff::Angles> angle_vec = {co2_angle1};
    eff::Conn connectivity =
        eff::Conn("co2_archive", bond_vec, angle_vec, std::vector<eff::Dihedral>{}, std::vector<eff::Inversion>{});
    std::unordered_map<std::string, eff::Conn> connectivity_map = {std::make_pair("co2_archive", connectivity)};

    // Now we create a system that will be the same as the one read
    bblock::System my_system;

    // Add monomers to the system
    size_t count = 0;
    for (size_t i = 0; i < n_monomers; i++) {
        std::vector<double> xyz(real_coords.begin() + 3 * count,
                                real_coords.begin() + 3 * count + 3 * n_atoms_vector[i]);
        std::vector<std::string> ats(atom_names.begin() + count, atom_names.begin() + count + n_atoms_vector[i]);
        std::string monid = monomer_names[i];
        my_system.AddMonomer(xyz, ats, monid);
        count += n_atoms_vector[i];
    }

    // Add the connectivity to the system
    my_system.SetConnectivity(connectivity_map);

    // Initialize the system to fill in the information
    my_system.Initialize();

    return my_system.OneBodyEnergy(false);


}




