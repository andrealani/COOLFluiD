#ifndef RANDOM_H
#define RANDOM_H

#include <cmath>

/**
 * Provides capabilities for computing random numbers.
 */
class Random
{
public:

    /**
     * Seeds the random number generator with the seed value.  If a see of 0 is
     * given, then the time will be used as a seed.
     */
    static void seed(unsigned int s = 0);

    /// Returns a uniformly distributed value between 0 and 1.
    static double uniform();

    /// Returns a uniformly distributed value between min and max.
    static double uniform(double min, double max) {
        return uniform()*(max - min) + min;
    }

    /// Returns a normally distributed number with mean 0 and deviation 1.
    static double normal();

    /// Returns a normally distributed number with mean and standard deviation.
    static double normal(double mean, double dev) {
        return normal()*dev + mean;
    }

    /// Returns a random number with a Voigt distribution given the Lorentz and
    /// Doppler half-widths at half-maximum (HWHM) and the center wavenumber.
    static double voigt(double gl, double gd, double sigc = 0.0)
    {
        static const double PI = 4.0 * std::atan(1.0);
        static const double RTLN2 = std::sqrt(std::log(2.0));
        return (gl*std::tan((uniform()-0.5)*PI) +
                gd/RTLN2*std::sqrt(-std::log(uniform()))*std::cos(2*PI*uniform()) + sigc);
    }

    /// Computes a uniformly distributed, random unit vector.
    template <typename Vector>
    static Vector unitVector();

private:

    /// Don't allow creation of Random objects.
    Random() { }

}; // class Random


template <typename Vector>
Vector Random::unitVector()
{
    Vector v;
    double sum = 0.0;
    for (int i = 0; i < v.size(); ++i) {
        v[i] = normal();
        sum += v[i]*v[i];
    }

    sum = 1.0/std::sqrt(sum);
    for (int i = 0; i < v.size(); ++i)
        v[i] *= sum;

    return v;
}


#endif // RANDOM_H
