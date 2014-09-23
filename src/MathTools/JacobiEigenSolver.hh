// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MathTools_JacobiEigenSolver_hh
#define COOLFluiD_MathTools_JacobiEigenSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "MatrixEigenSolver.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a eigen value solver for symmetric matrices
/// based on Jacobi algorithm
/// @author Jurek Majewski
class MathTools_API JacobiEigenSolver : public MatrixEigenSolver
{
public:

  /// Default constructor without arguments
  JacobiEigenSolver() : MatrixEigenSolver() {}

  /// Default destructor
  ~JacobiEigenSolver()  {}

  /// Finds an eigenvalues and eigenvectors of matrix a
  /// @param a    input matrix - it will be modified during calculation
  /// @param r    result: left eigenvectors matrix
  /// @param lam  result: vector of eigenvalues
  void eigenCalc( RealMatrix& a, RealMatrix& r, RealVector& lam );

private:
  void  RotateMatrix( RealMatrix& a,
    CFreal& g, CFreal& h, CFreal& s, CFreal& c, CFreal& tau,
    const CFuint& i, const CFuint& j, const CFuint& k, const CFuint& l );


}; // end of class JacobiEigenSolver

//////////////////////////////////////////////////////////////////////////////

  } // namespace MathTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MathTools_JacobiEigenSolver_hh
