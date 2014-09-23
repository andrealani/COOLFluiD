// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_EigenSolver_hh
#define COOLFluiD_Framework_EigenSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents an abstract eigensystem solver
  /// @author Andrea Lani
class Framework_API EigenSolver : public Common::OwnedObject {
public:

  typedef Environment::ConcreteProvider<EigenSolver,1> PROVIDER;
  typedef const CFuint ARG1;

  /// Default constructor without arguments
  EigenSolver(const CFuint n) :
    OwnedObject(),
    _n(n)
  {
  }

  /// Default destructor
  virtual ~EigenSolver()
  {
  }

  /// Computes the eigen values and right eigen vectors
  /// for a given non symmetric matrix with real coefficients
  /// @param A given matrix
  /// @param d array of the eigen values
  /// @param R matrix of right eigenvectors
  virtual void eigValuesVectorsRealNonSymm(const RealMatrix& A,
  				   RealVector& d,
  				   RealMatrix& R) = 0;

  /// Get the class name
  static std::string getClassName()
  {
    return "EigenSolver";
  }

protected: //data

  /// size of the matrix
  CFuint _n;

}; // end of class EigenSolver

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(EigenSolver) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_EigenSolver_hh
