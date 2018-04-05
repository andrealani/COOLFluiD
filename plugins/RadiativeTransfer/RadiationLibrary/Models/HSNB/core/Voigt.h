#ifndef COOLFluiD_RadiativeTransfer_VOIGT_H
#define COOLFluiD_RadiativeTransfer_VOIGT_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Constants.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LblSpectralGrid.h"
#include <cmath>

#include <iostream>


/**
 * Singleton class providing algorithms for computing the Voigt function for
 * line profiles.
 */
class Voigt
{
public:

    /**
     * Computes the Voigt function using the method of Humlicek.  The
     * implementation is converted from Fortran subroutine provedid in
     * J. Humlicek, J. Quant. Spectr. Rad. Transfer, v.27, pp.437-444, 1982.
     *
     * @param bl Lorentz half-width at half-maximum
     * @param bd Doppler half-width at half-maximum
     * @param ds Distance from the line center
     */
    static double humlicek(double bl, double bd, double ds);

    /**
     * Computes the Voigt function using the method of Drayson.  This
     * implementation has been converted to C++ from the EM2C Fortran version.
     *
     * @param bl Lorentz half-width at half-maximum
     * @param bd Doppler half-width at half-maximum
     * @param ds Distance from the line center
     */
    static double drayson(double bl, double bd, double ds);

    template <typename DataType, typename OP>
    static void drayson(
        const LblSpectralGrid& grid, DataType& data, double* const p_voigt,
        const OP& op);

    template<typename OP>
    static void drayson(
        const LblSpectralGrid& grid, double gamd, double gaml, double center,
        double* const p_voigt, const OP& op);

private:

    // Singleton
    Voigt() { }
};

//long r1 = 0;
//long r2 = 0;
//long r3a = 0;
//long r3b = 0;

double new_drayson(const double& x, const double& y);


/**
 * Computes the K(x,y) function using Drayson's method.  This implementation is
 * converted to C++ from the EM2C Fortran version.
 */
double drayson_helper(const double& x, const double& y);



/**
 * Computes the Voigt function over the SpectralGrid between indices start and
 * end, given the Doppler and Lorentz half-width at half-maximums.
 */
template <typename DataType, typename OP>
void Voigt::drayson(
    const LblSpectralGrid& grid, DataType& data, double* const p_voigt,
	const OP& op)
{
	drayson(grid, data.gamd, data.gaml, data.sigc, p_voigt, op);
}

/**
 * Computes the Voigt function over the SpectralGrid between indices start and
 * end, given the Doppler and Lorentz half-width at half-maximums.
 */
template <typename OP>
void Voigt::drayson(
    const LblSpectralGrid& grid, double gamd, double gaml, double center,
    double* const p_voigt, const OP& op)
{
    // Get the start and end indices to loop over [start,end)
    const double delta = std::sqrt(gamd*gamd+gaml*gaml)*100000.0;
    const int    start = grid.index(center-delta);
    const int    end   = grid.index(center+delta);

    const double sqrtln2 = std::sqrt(LN2);
    const double sqrtpi  = std::sqrt(PI);
    const double utogd   = sqrtln2/gamd;
    const double y = utogd*gaml;
    //const double y = gaml/gamd;
    double x;

    // Special case for y >= 11 (all x in Region IIIb)
    if (y > 11.0) {
        const double y2 = y*y;
        const double hh2 = 0.2820948;
        const double xx2 = 0.7071068;
        const double fac = utogd/sqrtpi*y;
        //const double fac = y/(sqrtpi*gamd);
        double u, v;
        //r3b += end-start;
        for (int i = start; i < end; ++i) {
            x = std::abs(utogd*(grid[i]-center));
            //x = std::abs((grid[i]-center)/gamd);
            u = x - xx2;
            v = x + xx2;
            op(p_voigt, fac*(hh2/(y2+u*u) + hh2/(y2+v*v)), i);
        }
    } else {
        for (int i = start; i < end; ++i) {
            x = std::abs(utogd*(grid[i]-center));
            //x = std::abs((grid[i]-center)/gamd);
            op(p_voigt, utogd/sqrtpi*new_drayson(x, y), i);
            //op(p_voigt, new_drayson(x,y)/(sqrtpi*gamd), i);
        }
    }

    //std::cout << "RI: " << r1 << ", RII: " << r2 << ", RIIIa: " << r3a << ", RIIIb: " << r3b << std::endl;
}


#endif
