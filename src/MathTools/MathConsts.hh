// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MathConsts_hh
#define COOLFluiD_MathTools_MathConsts_hh

//////////////////////////////////////////////////////////////////////////////

#include <limits>    // for std::numeric_limits

#include "Common/COOLFluiD.hh"
#include "Common/NonInstantiable.hh"
#include "MathTools/MathTools.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// Provides an a set of static functions for mathematical constants
/// @author Tiago Quintino
class MathConsts : public Common::NonInstantiable<MathConsts> {
public:
  /// Returns the maximum number representable with the chosen precision
  static CFreal CFintMax () { return std::numeric_limits<CFint>::max(); }
  /// Definition of the minimum number representable with the chosen precision.
  static CFreal CFintMin () { return std::numeric_limits<CFint>::min(); }
  /// Returns the maximum number representable with the chosen precision
  static CFuint CFuintMax () { return std::numeric_limits<CFuint>::max(); }
  /// Definition of the minimum number representable with the chosen precision.
  static CFuint CFuintMin () { return std::numeric_limits<CFuint>::min(); }
  /// Returns the maximum number representable with the chosen precision
  static CFreal CFrealMax () { return std::numeric_limits<CFreal>::max(); }
  /// Definition of the minimum number representable with the chosen precision.
  static CFreal CFrealMin () { return std::numeric_limits<CFreal>::min(); }
  /// Definition of the maximum difference recognazible between two numbers with
  /// the chosen precision. Usefull for comparisons to zero  with real numbers:
  /// @code std::abs(x) > MathTools::MathConsts::CFrealEps()  @endcode
  static CFreal CFrealEps () { return std::numeric_limits<CFreal>::epsilon(); }
  /// Definition of Infinity
  static CFreal CFrealInf () { return std::numeric_limits<CFreal>::infinity(); }
  /// Definition of the Pi constant.
  static CFreal CFrealPi  () { return M_PI; }
  /// Definition of the imaginary constant i = sqrt(-1)
  static CFcomplex CFcomplexI () { return CFcomplex(0.0,1.0); }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_MathConsts_hh
