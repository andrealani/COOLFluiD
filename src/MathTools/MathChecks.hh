// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MathChecks_hh
#define COOLFluiD_MathTools_MathChecks_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonInstantiable.hh"

#include "MathTools/MathConsts.hh"

/*
#if (defined(macintosh) || defined(__APPLE__) || defined(__APPLE_CC__))
#  define LDBL_MANT_DIG -1 // avoid anoying warning in boost 1.42 on APPLE
#endif
*/

#ifndef CF_HAVE_CUDA
#include "boost/math/special_functions/fpclassify.hpp"
#endif

/// Define the cfFinite and cfIsNaN macros
#if defined(_MSC_VER)     /* Microsoft Visual C++ */
    #include <float.h>   // for _isnan
    #define cfFinite(n) _finite(n)
    #define cfIsNaN(n)  _isnan(n)
#elif defined(__GNUC__) /* GNU C++ */
    #define cfFinite(n) std::isfinite(n)
    #define cfIsNaN(n)  std::isnan(n)
#else // other compilers
    #define cfFinite(n) 1 // don't know
    #define cfIsNaN(n) ((n) != (n))
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// Provides an a set of static functions for checking CFreal numbers.
/// This class is not instantiable.
/// @author Tiago Quintino
class MathChecks : public Common::NonInstantiable<MathChecks> {
public:

  /// Function to check if two CFreal numbers are equal.
  /// @param x
  /// @param y
  /// @return true if equal or almost equal within the accepted error fuzz.
  static bool isEqual(const CFreal& x,  const CFreal& y)
  {
    return isEqualWithError(x,y,MathTools::MathConsts::CFrealMin());
  }

  /// Function to check if two CFreal numbers are not equal.
  /// @param x
  /// @param y
  /// @return true if not equal or almost unequal within the accepted error fuzz.
  static bool isNotEqual(const CFreal& x,  const CFreal& y)
  {
    return isNotEqualWithError(x,y,MathTools::MathConsts::CFrealMin());
  }

  /// Function to check if two CFreal numbers are equal.
  /// @param x
  /// @param y
  /// @param fuzz
  /// @return true if equal or almost equal within the accepted error fuzz.
  static bool isEqualWithError(const CFreal& x, const CFreal& y, const CFreal& fuzz)
  {
    // see Knuth section 4.2.2 pages 217-218
    return std::abs(x - y) <= fuzz * std::abs(x);
  }

  /// Function to check if two CFreal numbers are not equal.
  /// @param x
  /// @param y
  /// @param fuzz
  /// @return true if not equal or almost unequal within the accepted error fuzz.
  static bool isNotEqualWithError(const CFreal& x, const CFreal& y, const CFreal& fuzz)
  {
    return !isEqualWithError(x,y,fuzz);
  }

  /// Function to check if a CFreal number is finite number. This means is not a NaN neither a INF
  /// @param x
  /// @return true if x is finite
 //  static bool isFinite(const CFreal& x)
//   {
//     return (boost::math::isfinite)(x);
//   }

#ifndef CF_HAVE_CUDA
  /// Function to check if a CFreal number is either minus or plus INF
  /// @param x
  /// @return true if x is finite
  static bool isInf(const CFreal& x)
  {
    return (boost::math::isinf)(x);
  }

  /// Function to check if a CFreal number is a NaN (Not a Number)
  /// @param x
  /// @return true if x is a NaN
  static bool isNaN(const CFreal& x)
  {
    return (boost::math::isnan)(x);
  }
#else
 /// Function to check if a CFreal number is either minus or plus INF
  /// @param x
  /// @return true if x is finite
  static bool isInf(const CFreal& x) {return false;}

  /// Function to check if a CFreal number is a NaN (Not a Number)
  /// @param x
  /// @return true if x is a NaN
  static bool isNaN(const CFreal& x) {return false;}
#endif

  /// Function to check if a CFreal number is zero or very close.
  /// @param x
  /// @return true if equal to zero or almost equal within the accepted error fuzz.
  static bool isZeroWithError(const CFreal& x, const CFreal& fuzz)
  {
    return std::abs(x) <= fuzz;
  }


  /// Function to check if a CFreal number is not zero or very close.
  /// @param x
  /// @return true if not equal to zero or not almost equal within the accepted error fuzz.
  static bool isNotZeroWithError(const CFreal& x, const CFreal& fuzz)
  {
    return !isZeroWithError(x,fuzz);
  }

  /// Function to check if a CFreal number is zero or very close.
  /// @param x
  /// @return true if equal to zero or almost equal within the accepted error fuzz.
  static bool isZero(const CFreal& x)
  {
    return isZeroWithError(x,MathTools::MathConsts::CFrealMin());
  }

  /// Function to check if a CFreal number is not zero or very close.
  /// @param x
  /// @return true if not equal to zero or not almost equal within the accepted error fuzz.
  static bool isNotZero(const CFreal& x)
  {
    return isNotZeroWithError(x,MathTools::MathConsts::CFrealMin());
  }

  /// Sign function returning a real
  static CFreal sign(const CFreal& value)
  {
    return (value < 0.0) ? -1.0 : 1.0;
  }

  /// Checks is real is positive.
  /// Kind of a sign function returning a bool.
  static bool isPositive(const CFreal& value)
  {
    return (value < 0.0) ? false : true;
  }

  /// Checks is real is negative.
  /// Kind of a sign function returning a bool.
  static bool isNegative(const CFreal& value)
  {
    return !isPositive(value);
  }

}; // end class MathChecks

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_MathChecks_hh
