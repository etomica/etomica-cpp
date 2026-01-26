//
// Created by skwab on 12/6/2025.
//

#ifndef MYPROJECT_TRAPPEPARAMS_H
#define MYPROJECT_TRAPPEPARAMS_H
#include "P2PotentialGroupBuilder.h"
#include "potential-master.h"
#include "species.h"


class TraPPEParams {
public:
    enum ChemForm {
        CO2, propane
        };
    TraPPEParams(ChemForm chemForm);
    ChemForm chemForm;
    vector<double> sigma;
    vector<double> epsilon;
    vector<double> charge;
    vector<double> k_theta;
    vector<double> theta_eq;
    vector<double> k_b;
    vector<double> r_eq;
    vector<vector<double>> a;
    vector<vector<int *>> pairs;
    vector<int *> triplets;
    vector<vector<int>> quads;
    vector<vector<int>> bonding;
    SpeciesSimple species;
    SpeciesList speciesList;
    string diagram;
    // Box box;
    bool isFlex = false;
    bool polar = false;
private:
    SpeciesSimple makeSpecies(ChemForm chemForm);

};


#endif //MYPROJECT_TRAPPEPARAMS_H