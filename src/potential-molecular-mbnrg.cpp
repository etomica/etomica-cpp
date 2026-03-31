#include "potential-molecular.h"
#include "bblock/system.h"
#include "unit.h"
#include "vector.h"

PotentialMolecularMBnrg::PotentialMolecularMBnrg() : pot("co2_archive", "co2_archive"), grad(9, 0), chg(6, 0), chggrad(6, 0), pol(6, 0), polfac(6, 0)
{
    xyz1 = (double*)malloc(sizeof(double) * 3 * 3);
    xyz2 = (double*)malloc(sizeof(double) * 3 * 3);
    vector<double> sys_xyz(6*3, 0.0);
    nlohmann::json monomers_j_{};
    vector<size_t> first_index = {0, 0};
    systools::SetCharges(sys_xyz, chg, "co2_archive", 2, 3, 0, chggrad, monomers_j_);
    systools::SetPolfac(polfac, "co2_archive", 2, 3, 0, monomers_j_);
    systools::SetPol(pol, "co2_archive", 2, 3, 0, monomers_j_);
    vector<size_t> sites = {3, 3};
    vector<string> mon_id = {"co2_archive", "co2_archive"};
    vector<size_t> num_atoms = {3, 3};
    vector<pair<string, size_t>> mon_type_count = {{"co2_archive", 2}};
    vector<size_t> islocal = {1, 1};
    vector<int> tags = {0, 1};
    vector<double> edge_vectors;
    vector<double> c6_lr(6, 0);

    electrostatics.Initialize(chg, chggrad, polfac, pol, sys_xyz, mon_id, sites, first_index, mon_type_count,
                               islocal, tags, false, 1E-16, 100, "iter", edge_vectors);
    systools::SetC6LongRange(c6_lr, "co2_archive", 2, 3, 0, monomers_j_);


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
    double minDist = 1000;
    for (int j=3; j<6; j++)
    {
        double dr2[3];
        for (int l=0;l<3;l++)
        {
            for (int k=0; k<3; k++) dr2[k] = xyz2[3*(j-3)+k]-xyz1[3*(l)+k];
            double r2_2 = Vector::squared(dr2);
            minDist = min(minDist, r2_2);
        }
    }
    if (minDist < 1) return std::numeric_limits<double>::infinity();


    double energy2b = pot.eval(xyz1, xyz2, 1);
    double total_energy = energy2b;
    if (energy2b < 100000)
    {
        vector<double> edge_vectors;
        vector<pair<string, string>> ignore_disp;
        dispersion.SetNewParameters(xyz, ignore_disp, false, std::numeric_limits<double>::infinity(), edge_vectors);
        double disp_energy = dispersion.GetDispersion(grad, nullptr, false);
        electrostatics.SetNewParameters(xyz, chg, chggrad, pol, polfac, "iter", false, edge_vectors, std::numeric_limits<double>::infinity());
        double elec_energy = electrostatics.GetElectrostatics(grad, nullptr, false);
        total_energy += elec_energy + disp_energy;
        double badEnergy = -10;
        double highEnergy = 400;
        if (energy2b < badEnergy)
        {
            for (int x = 0; x <2 ; x++)
            {
                double* p = box.getAtomPosition(0);
                double* q = box.getAtomPosition(3);
                double dr[3];
                for (int k=0; k<3; k++) dr[k] = q[k]-p[k];
                double r2 = sqrt(Vector::squared(dr));
                for (int k=0; k<3; k++) dr[k]/=r2;
                double minDistHighEnergy = 100;
                double maxDistHighEnergy = 0;
                double step = 0.05;
                for (int i=0; i<100&&r2>0.01; i++)
                {
                    double minDist = 1000;
                    for (int j=3; j<6; j++)
                    {
                        for (int k=0; k<3; k++) xyz2[3*(j-3)+k]+=step*dr[k];
                        if (isnan(xyz2[3*(j-3)]))
                        {
                            printf("%d %f %d\n",  j, box.getAtomPosition(j)[0], i);
                            exit(0);
                        }
                        double dr2[3];
                        for (int l=0;l<3;l++)
                        {
                            for (int k=0; k<3; k++) dr2[k] = xyz2[3*(j-3)+k]-xyz1[3*(l)+k];
                            double r2_2 = Vector::squared(dr2);
                            minDist = min(minDist, r2_2);
                        }
                    }
                    double u2 = pot.eval(xyz1, xyz2, 1);
                    double dist = sqrt(minDist);

                    if (x==1) printf("dist = %f, U = %f\n", dist, u2);
                    if (x==1&&u2<highEnergy&&dist<0.9)
                    {
                        printf("6\n");
                        printf("%f \n", total_energy);
                        printf("C %f %f %f\n", xyz1[0], xyz1[1], xyz1[2]);
                        printf("O %f %f %f\n", xyz1[3], xyz1[4], xyz1[5]);
                        printf("O %f %f %f\n", xyz1[6], xyz1[7], xyz1[8]);

                        printf("C %f %f %f\n", xyz2[0], xyz2[1], xyz2[2]);
                        printf("O %f %f %f\n", xyz2[3], xyz2[4], xyz2[5]);
                        printf("O %f %f %f\n", xyz2[6], xyz2[7], xyz2[8]);

                        break;
                    }
                    if (u2 > highEnergy && dist < minDistHighEnergy)
                    {
                        minDistHighEnergy = dist;
                        maxDistHighEnergy = dist;
                    }
                    else if (u2 > highEnergy && dist > maxDistHighEnergy)
                    {
                        maxDistHighEnergy = dist;

                    }
                    else if (u2 < badEnergy)
                    {
                        minDistHighEnergy = 100;
                    }
                }
                bool bad = false;
                if (r2 > 0.01&&x==0)
                {
                    if (minDistHighEnergy == 100)
                    {
                        printf("Never encountered high energy");
                        exit(0);
                    }
                    if (maxDistHighEnergy > 100*step)
                    {
                        printf("Never encountered low energy");
                        exit(0);
                    }
                    bool printme = maxDistHighEnergy < minMax || minDistHighEnergy > maxMin;
                    minMax = min(minMax, maxDistHighEnergy);
                    maxMin = max(maxMin, minDistHighEnergy);
                    if (printme) printf("minMax = %g, maxMin = %g\n", minMax, maxMin);
                    if (maxDistHighEnergy < 0.9)
                    {
                        for (int j = firstAtom; j <= lastAtom; j++)
                        {
                            double* p = box.getAtomPosition(j);
                            xyz2[3*(j-firstAtom)] = p[0];
                            xyz2[3*(j-firstAtom)+1] = p[1];
                            xyz2[3*(j-firstAtom)+2] = p[2];
                        }
                        bad = true;
                    }

                }
                if (!bad) break;
            }

            // printf("6\n");
            // printf("%f \n", total_energy);
            // printf("C %f %f %f\n", xyz1[0], xyz1[1], xyz1[2]);
            // printf("O %f %f %f\n", xyz1[3], xyz1[4], xyz1[5]);
            // printf("O %f %f %f\n", xyz1[6], xyz1[7], xyz1[8]);
            //
            // printf("C %f %f %f\n", xyz2[0], xyz2[1], xyz2[2]);
            // printf("0 %f %f %f\n", xyz1[3], xyz1[4], xyz1[5]);
            // printf("O %f %f %f\n", xyz1[6], xyz1[7], xyz1[8]);
            total_energy = 1000000;


        }
    }
    return Energy::toSim(total_energy);
}