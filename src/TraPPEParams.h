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
    TraPPEParams();

    void buildPotentials();

    int numAtomTypes;
    vector<double> sigma;
    vector<double> epsilon;
    vector<double> charge;
    vector<double> k_theta;
    vector<double> theta_eq;
    vector<vector<double>> a;
    vector<vector<int>> triplets;
    vector<vector<int>> quads;
    vector<vector<int>> bonding;
    SpeciesSimple species;
    SpeciesList speciesList;
    // Box box;
    enum chemForm {
        CO2
        };

    chemForm chemForm;
    // struct params {
    //     int numAtomTypes;
    //     vector<double> sigma;
    //     vector<double> epsilon;
    //     vector<double> charge;
    // };
    // params params;
    // PotentialMasterVirial potentialGroup;
    bool isFlex = false;
    bool polar = false;


};


#endif //MYPROJECT_TRAPPEPARAMS_H