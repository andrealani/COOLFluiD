// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MatrixIntersect_hh
#define COOLFluiD_MathTools_MatrixIntersect_hh

//////////////////////////////////////////////////////////////////////////////

#include "RealMatrix.hh"
#include "RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class is used to find an intersection of two SPD matrices
   *
   * @author Jurek Majewski
   */
class MathTools_API MatrixIntersect
{
public:

  /**
   * Default constructor without arguments
   */
  MatrixIntersect()  {}

  /**
   * Default destructor
   */
  virtual ~MatrixIntersect()	{}

  /**
   * Factory method to create an inverter suitable for the
   * given size
   */
  static MatrixIntersect* create( const CFuint& size);

  /**
   * Invert the given matrix a and put the result in x
   * @param a    first matrix
   * @param b    second matrix
   * @param res  result: final matrix
  */
  virtual void intersectCalc( const RealMatrix& a, const RealMatrix& b, RealMatrix& res) = 0;

}; // end of class MatrixIntersect

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_MatrixIntersect_hh
