// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_LinearAdvSys3DPrim_hh
#define COOLFluiD_Physics_LinearAdvSys_LinearAdvSys3DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinearAdvSys/LinearAdvSys3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Linearized Euler physical model 3D for conservative
 * variables
 *
 * @author Tiago Quintino
 * @author Tomas Kopacek
 */
class LinearAdvSys3DPrim : public LinearAdvSys3DVarSet {
public: // classes

  /**
   * Constructor
   * @see LinearAdvSys3D
   */
  LinearAdvSys3DPrim(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~LinearAdvSys3DPrim();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

   /**
   * Set the jacobian matrix
   */
  virtual void computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob);

  /**
   * Split the jacobian
   */
  void splitJacobian(RealMatrix& jacobPlus,
      RealMatrix& jacobMin,
      RealVector& eValues,
      const RealVector& normal);

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  void computeEigenValuesVectors(RealMatrix& rightEv,
           RealMatrix& leftEv,
           RealVector& eValues,
           const RealVector& normal);

  /**
   * Set the first right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  void setEigenVect1(RealVector& r1,
         Framework::State& state,
         const RealVector& normal);
  /**
   * Set the second right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  void setEigenVect2(RealVector& r2,
         Framework::State& state,
         const RealVector& normal);
  /**
   * Set the third right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  void setEigenVect3(RealVector& r3,
         Framework::State& state,
         const RealVector& normal);

  /**
   * Set the fourth right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  void setEigenVect4(RealVector& r4,
         Framework::State& state,
         const RealVector& normal);

  /**
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const;

  /**
   * Give dimensional values to the adimensional state variables
   */
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  void setAdimensionalValues(const Framework::State& state,
                             RealVector& result);

  /**
   * Compute the perturbed states data
   */
  void computePerturbedStatesData
  (const std::vector<Framework::State*>& states,
   const CFuint nbStatesInVec,
   const CFuint iVar);

  /**
   * Set the PhysicalData corresponding to the given State
   * @see LinearAdvSysPhysicalModel
   */
  /// @todo broken after release 2009.3 (added "const")
  void computePhysicalData(const Framework::State& state,
  		   RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see LinearAdvSysPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data,
			   Framework::State& state);
  
protected:
  
  /**
   * Set the constant part (independent form the solution) of the
   * jacobians
   */
  void setConstJacob();


private: // data

  /// temporary matrix of right eigenvalues
  RealMatrix                       _rightEv;

  /// temporary matrix of left eigenvalues
  RealMatrix                       _leftEv;

}; // end of class LinearAdvSys3DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdvSys_LinearAdvSys3DPrim_hh
