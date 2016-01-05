// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_FunctionParser_hh
#define COOLFluiD_MathTools_FunctionParser_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "MathTools/FParser/fparser.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a simple interface to the Warp's function parser 
/// @see website http://warp.povusers.org/FunctionParser/
/// @author Andrea Lani

class MathTools_API FunctionParser_f: public FunctionParserBase<float> {};
class MathTools_API FunctionParser_ld: public FunctionParserBase<long double> {};
class MathTools_API FunctionParser_li: public FunctionParserBase<long> {};

#include <complex>
class MathTools_API FunctionParser_cd: public FunctionParserBase<std::complex<double> > {};
class MathTools_API FunctionParser_cf: public FunctionParserBase<std::complex<float> > {};
class MathTools_API FunctionParser_cld: public FunctionParserBase<std::complex<long double> > {};

//----------------------------------------------------------------------------//

/// New functions to be added to the default ones 

/// radius in 2D
static double Radius2D(const double* x)
{
  return std::sqrt(x[0]*x[0] + x[1]*x[1]);
} 

/// radius in 3D
static double Radius3D(const double* x)
{
  return std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

/// Random function with positive numbers from [0,1]
static double Rand(const double*)
{
    return drand48();
} 
    
//----------------------------------------------------------------------------//

// This class derives from the FunctionParser defined in the library from 
/// @see website http://warp.povusers.org/FunctionParser/
// DerivedParser adds some constants and functions to the set already available 
// in the original library

// @author Andrea Lani

class MathTools_API FunctionParser: public FunctionParserBase<double> {
public:
  /// default constructor
  FunctionParser() : FunctionParserBase<double>()
  {
    AddConstant("pi", 3.14159265358979323846);
    AddConstant("e", 2.71828182845904523536);
    AddConstant("Rair", 287.046);
    AddFunction("R2", Radius2D, 2);
    AddFunction("R3", Radius3D, 3);
    AddFunction("rand", Rand, 0);
  }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_FunctionParser_hh
