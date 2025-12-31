#include "unit.h"
#include "constants.h"





double Electron::convFactor = 1.0 / std::sqrt(4.0 * M_PI * Constants::EPSILON_0);

double Degree::convFactor = M_PI/180;