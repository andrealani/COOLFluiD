#ifndef COOLFluiD_RadiativeTransfer_CONSTANTS_H
#define COOLFluiD_RadiativeTransfer_CONSTANTS_H

#include <cmath>
#include <limits>

// Pi
const double PI     = 4.0 * std::atan(1.0);
const double TWOPI  = 2.0 * PI;
const double FOURPI = 4.0 * PI;

// Universal constants (all units in SI)
const double NA     = 6.0221415E23;    // Avagadro's number (molecule/mol)
const double KB     = 1.380658E-23;   // Boltzmann's constant (J/molecule-K)
const double RU     = NA * KB;         // Universal Gas constant (J/mole-K)
const double HP     = 6.6260755E-34;    // Planck's constant (J-s)
const double MU0    = PI * 4.0E-7;     // Magnetic constant (H/m)
const double C0     = 299792458.0;     // Speed of light in vacuum (m/s)
const double EPS0   = 1.0/(MU0*C0*C0); // Vacuum permittivity (F/m)
const double QE     = 1.6021773349E-19; // Elementary positive charge (C)
const double ONEATM = 101325.0;        // 1 atm in Pa
const double ONEBAR = 100000.0;        // 1 bar in Pa
const double EIONH  = 2.1798741E-18;   // Ionization energy of hydrogen (J/atom)
const double ME     = 9.109389754E-31;  // Electron mass (kg)

// Common values
const double SQRT2  = std::sqrt(2.0);  // square root of 2
const double LN2    = std::log(2.0);   // natural log of 2

const double MAXLOG = std::log(std::numeric_limits<double>::max());

const double F1TIMESSIGSQUARED = 0.25e-4*QE*QE/(C0*C0*ME*EPS0);


#endif // CONSTANTS_H
