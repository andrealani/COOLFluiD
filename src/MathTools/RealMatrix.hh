// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_RealMatrix_hh
#define COOLFluiD_RealMatrix_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CXX_EXPLICIT_TEMPLATES
  // explicit template instantiation
  MathTools_TEMPLATE template class MathTools_API MathTools::CFMat<CFreal>;
  MathTools_TEMPLATE template class MathTools_API MathTools::CFMatSlice<CFreal>;
#endif

/// RealMatrix is a CFMat templatized with a CFreal
typedef MathTools::CFMat<CFreal> RealMatrix;

/// RealSliceMatrix is a CFSliceMatrix templatized with a CFreal
typedef MathTools::CFMatSlice<CFreal> RealSliceMatrix;

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RealMatrix_hh
