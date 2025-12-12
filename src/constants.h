//
// Created by skwab on 12/6/2025.
//

#ifndef MYPROJECT_CONSTANTS_H
#define MYPROJECT_CONSTANTS_H
#include <cmath>

class Constants {
    public:
    static constexpr double AVOGADRO = 6.02214076e23;
    static constexpr double EPSILON_0 =
        8.8541878128e-12 * (1.0 / (1.602176634e-19 * 1.602176634e-19)) *
        1e24 * 1e-3 / AVOGADRO * 1e-30;

    static constexpr double BOLTZMANN_K =
        1.380649e-23 * 1000 * AVOGADRO * 1e20 * 1e-24;
};

#endif //MYPROJECT_CONSTANTS_H