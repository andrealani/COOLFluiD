// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_MatrixEigenSolver_hh
#define COOLFluiD_MathTools_MatrixEigenSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "RealMatrix.hh"
#include "RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a matrix eigen value problem solver
  /// @author Jurek Majewski
class MathTools_API MatrixEigenSolver
{
public:

  /// Default constructor without arguments
  MatrixEigenSolver()  {}

  /// Default destructor
  virtual ~MatrixEigenSolver()  {}

  /// Factory method to create an inverter suitable for the
  /// given size
  static MatrixEigenSolver* create(const CFuint& size, const bool& isSymmetric);

  /// Invert the given matrix a and put the result in x
  /// @param a    matrix to invert
  /// @param r    result: left eigenvectors matrix
  /// @param lam  result: vector of eigenvalues
  virtual void eigenCalc( RealMatrix& a, RealMatrix& r, RealVector& lam) = 0;

}; // end of class MatrixEigenSolver

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_MatrixEigenSolver_hh
