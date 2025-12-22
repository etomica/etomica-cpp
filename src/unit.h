//
// Created by skwab on 12/6/2025.
//

#ifndef MYPROJECT_UNIT_H
#define MYPROJECT_UNIT_H
#include <math.h>

#include "constants.h"

class Kelvin{
    public:
    static constexpr double convFactor = Constants::BOLTZMANN_K;

    static double fromSim(double x) { return x / convFactor; }
    static double toSim(double x)   { return x * convFactor; }
};

class Electron{
    public:
    static double convFactor;


    static double fromSim(double x) { return x / convFactor; }
    static double toSim(double x)   { return x * convFactor; }
};


#endif //MYPROJECT_UNIT_H