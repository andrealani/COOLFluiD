// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_RealVector_hh
#define COOLFluiD_RealVector_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFVec.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CXX_EXPLICIT_TEMPLATES
  // explicit template instantiation
  MathTools_TEMPLATE template class MathTools_API MathTools::CFVec<CFreal>;
  MathTools_TEMPLATE template class MathTools_API MathTools::CFVecSlice<CFreal>;
#endif

/// RealVector is a CFVec templatized with a CFreal
typedef MathTools::CFVec<CFreal> RealVector;

/// RealSliceVector is a CFSliceVector templatized with a CFreal
typedef MathTools::CFVecSlice<CFreal> RealSliceVector;

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RealVector_hh
