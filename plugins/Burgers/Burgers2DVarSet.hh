// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_Burgers_Burgers2DVarSet_hh
#define COOLFluiD_Physics_Burgers_Burgers2DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Burgers physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 *
 */
class Burgers2DVarSet : public Framework::ConvectiveVarSet {
public: // classes

  /**
   * Constructor
   * @see Burgers2D
   */
  Burgers2DVarSet(Common::SafePtr<Framework::BaseTerm> term) :
    Framework::ConvectiveVarSet(term)
  {
  }

  /**
   * Default destructor
   */
  virtual ~Burgers2DVarSet()
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
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal) = 0;

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  /// @todo broken after release 2009.3 (added "const")
  virtual void computePhysicalData(const Framework::State& state, RealVector& data) = 0;

  /**
   * Set the PhysicalData corresponding to the given State
   */
  virtual void computeStateFromPhysicalData(const RealVector& data, Framework::State& state) = 0;

protected:

  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);

  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);

  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);

}; // end of class Burgers2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace Burgers

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Burgers_Burgers2DVarSet_hh
