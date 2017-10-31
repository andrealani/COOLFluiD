#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>

#include "Random.h"

//==============================================================================

void Random::seed(unsigned int s) {
    srand(s == 0 ? time(NULL) : s);
}

//==============================================================================

double Random::uniform() {
    return (rand() * (1.0 / RAND_MAX));
}

//==============================================================================

double Random::normal()
{
    // Using the Box-Muller tranform.
    // https://en.wikipedia.org/wiki/Box-Muller_transform
    const double epsilon = std::numeric_limits<double>::min();
    const double two_pi = 2.0*3.14159265358979323846;

    static double z0, z1;
    static bool generate = false;
    generate = !generate;

    if (!generate)
       return z1;

    double u1, u2;
    do {
       u1 = rand() * (1.0 / RAND_MAX);
    } while ( u1 <= epsilon );
    u2 = rand() * (1.0 / RAND_MAX);

    u1 = sqrt(-2.0 * log(u1));
    u2 *= two_pi;
    z0 = u1 * cos(u2);
    z1 = u1 * sin(u2);
    return z0;
}

//==============================================================================
