// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullLinearSystemSolver_hh
#define COOLFluiD_Framework_NullLinearSystemSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinearSystemSolver.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullLinearSystemSolver.
/// @author Andrea Lani
class Framework_API NullLinearSystemSolver : public LinearSystemSolver {
public:

  /// Default constructor without arguments
  NullLinearSystemSolver(const std::string& name);

  /// Default destructor
  ~NullLinearSystemSolver();

  /// Solve the linear system
  void solveSys();

  /// Prints the Linear System to a file.
  void printToFile(const std::string prefix, const std::string suffix);

  /// Create a block accumulator with chosen internal storage
  /// @return a newly created block accumulator
  /// @post the block has to be deleted outside
  BlockAccumulator* createBlockAccumulator(const CFuint nbRows,
					   const CFuint nbCols,
					   const CFuint subBlockSize,
					   CFreal* ptr) const;
  
  /// Get the LSS system matrix
  virtual Common::SafePtr<LSSMatrix> getMatrix() const;

  /// Get the LSS solution vector
  virtual Common::SafePtr<LSSVector> getSolVector() const;

  /// Gets the LSS right hand side vector
  virtual Common::SafePtr<LSSVector> getRhsVector() const;

  /// Checks if this object is a Null object.
  /// Since this is NullSpaceMethod
  /// @return true
  virtual bool isNull() const { return true; }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Solve the linear system
  /// @see LinearSystemSolver::solveSys()
  virtual void solveSysImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

}; // class NullLinearSystemSolver

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullLinearSystemSolver_hh
