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

class MathTools_API FunctionParser: public FunctionParserBase<double> {};
class MathTools_API FunctionParser_f: public FunctionParserBase<float> {};
class MathTools_API FunctionParser_ld: public FunctionParserBase<long double> {};
class MathTools_API FunctionParser_li: public FunctionParserBase<long> {};

#include <complex>
class MathTools_API FunctionParser_cd: public FunctionParserBase<std::complex<double> > {};
class MathTools_API FunctionParser_cf: public FunctionParserBase<std::complex<float> > {};
class MathTools_API FunctionParser_cld: public FunctionParserBase<std::complex<long double> > {};

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_FunctionParser_hh
