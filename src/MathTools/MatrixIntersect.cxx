// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "IntersectSolver.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

MatrixIntersect* MatrixIntersect::create( const CFuint& size )
{
  if ( size > 10 )
    std::cout << "MatrixIntersect::create() : size of matrix is greater then 10. a better algorithm should be used\n" << std::endl;

  return new IntersectSolver();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
