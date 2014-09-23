// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MatrixInverter_hh
#define COOLFluiD_MathTools_MatrixInverter_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a matrix inverter
/// @author Andrea Lani
/// @author Tiago Quintino
class MathTools_API MatrixInverter {
public:

  /// Default constructor without arguments
  MatrixInverter() {}

  /// Default destructor
  virtual ~MatrixInverter() {}

  /// Factory method to create an inverter suitable for the given size
  static MatrixInverter* create (const CFuint& size, const bool& isDiagonal);
  
  /// Factory method to create an inverter suitable for non square matrices (pseudo-inverse)
  static MatrixInverter* create (const CFuint& nbRows, const CFuint& nbCols, const bool& isDiagonal);

  /// Invert the given matrix a and put the result in x
  /// @param a  matrix to invert
  /// @param x  result of the matrix inversion
  virtual void invert(const RealMatrix& a, RealMatrix& x) = 0;

}; // end of class MatrixInverter

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_MatrixInverter_hh
