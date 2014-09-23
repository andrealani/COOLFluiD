// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_NonLinearAdv_NonLinearAdv2DPrim_hh
#define COOLFluiD_Physics_NonLinearAdv_NonLinearAdv2DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "NonLinearAdv2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {
      class NonLinearAdvTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NonLinearAdv physical model 2D for conservative
 * variables
 *
 * @author Nadege Villedieu
 *
 */
class NonLinearAdv2DPrim : public NonLinearAdv2DVarSet {
public:
  /**
   * Constructor
   * @see NonLinearAdv2D
   */
  NonLinearAdv2DPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~NonLinearAdv2DPrim();

  /**
   * Gets the average advection vector the specified variable
   * @param vec the vector to change
   * @param iVar the variable corresponding to the specified vector
   */
  virtual void getAverageAdvectionVector(RealVector& vec, const CFuint iVar) const;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Set the scalar part of the jacobian
   */
  virtual void computeScalarJacobian(const RealVector& normal,
          RealVector& jacob);

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void setEigenValuesVectors(RealMatrix& rightEv,
                             RealMatrix& leftEv,
                             RealVector& eValues,
                             const RealVector& normal);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  /// @todo broken after release 2009.3 (added "const")
  void computePhysicalData(const Framework::State& state,
			   RealVector& data);
 
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  { 
    velIDs.resize(1); velIDs[XX] = 1;
  }
  
protected:

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const Framework::State& vars,
			   const RealVector& normals);

  /**
   * Compute the physical convective flux
   */
  virtual void computeFlux(const Framework::State& vars);

  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);
  
  virtual void computeStateFlux(const RealVector& pdata);
  
}; // end of class NonLinearAdv2DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NonLinearAdv_NonLinearAdv2DPrim_hh
