// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_NonLinearAdv_NonLinearAdv2DVarSet_hh
#define COOLFluiD_Physics_NonLinearAdv_NonLinearAdv2DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "NonLinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NonLinearAdv physical model 2D for conservative
 * variables
 *
 * @author Nadege Villedieu
 *
 */
class NonLinearAdv2DVarSet : public Framework::ConvectiveVarSet {

public: // classes

  /**
   * Constructor
   * @see NonLinearAdv2D
   */
  NonLinearAdv2DVarSet(Common::SafePtr<Framework::BaseTerm> term) :
    Framework::ConvectiveVarSet(term),
    _model(term.d_castTo<NonLinearAdvTerm>())
  {
  }

  /**
   * Default destructor
   */
  virtual ~NonLinearAdv2DVarSet()
  {
  }

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();
  }

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const
  {
    return 1;
  }

  /**
   * Set the jacobians
   */
  virtual void computeJacobians() = 0;

  /**
   * Set the scalar part of the jacobian
   */
  virtual void computeScalarJacobian(const RealVector& normal,
                              RealVector& jacob) = 0;

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void setEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal) = 0;

  /**
   * Get the model
   */
  Common::SafePtr<NonLinearAdvTerm> getModel() const
  {
    return _model;
  }

protected:

  /**
   * Set the vector of the eigenValues
   */
  virtual void setEigenValuesImpl(Framework::State& state,
			  const RealVector& normal, RealVector& eValues);

  /**
   * Get the maximum eigen value
   */
  virtual CFreal getMaxEigenValueImpl(Framework::State& state, const RealVector& normal);

  /**
   * Get the maximum absolute eigenvalue
   */
  virtual CFreal getMaxAbsEigenValueImpl(Framework::State& state, const RealVector& normal);

private:

  /// acquaintance of the model
  Common::SafePtr<NonLinearAdvTerm> _model;

}; // end of class NonLinearAdv2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NonLinearAdv_NonLinearAdv2DVarSet_hh
