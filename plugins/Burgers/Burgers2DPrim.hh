// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_Burgers_Burgers2DPrim_hh
#define COOLFluiD_Physics_Burgers_Burgers2DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Burgers2DVarSet.hh"
#include "BurgersPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

      class BurgersTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Burgers physical model 2D for primitive
 * variables
 *
 * @author Andrea Lani
 *
 */
class Burgers2DPrim : public Burgers2DVarSet {
public:

  /**
   * Constructor
   * @see Burgers2D
   */
  Burgers2DPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Burgers2DPrim();

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Set the scalar part of the jacobian
   */
  void computeScalarJacobian(const RealVector& normal, RealVector& jacob);

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  void computeEigenValuesVectors(RealMatrix& rightEv,
			     RealMatrix& leftEv,
			     RealVector& eValues,
			     const RealVector& normal);
  /**
   * Get the model
   */
  Common::SafePtr<BurgersTerm> getModel() const
  {
    return _model;
  }

  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  void computePhysicalData(const Framework::State& state, RealVector& data);

  /**
   * Set the PhysicalData corresponding to the given State
   */
  void computeStateFromPhysicalData(const RealVector& data, Framework::State& state)
  {
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(1); velIDs[XX] = 1;
  }
  
protected:

  /**
   * Compute the convective flux
   */
  virtual void computeFlux(const RealVector& vars,
			   const RealVector& normals);
  
  virtual void computeFlux(const RealVector& vars);

  /**
   * Compute the physical convective flux
   */
  virtual void computeStateFlux(const RealVector& vars);

private:

  /// acquaintance of the model
  Common::SafePtr<BurgersTerm> _model;

}; // end of class Burgers2DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Burgers

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Burgers_Burgers2DPrim_hh
