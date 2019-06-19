#ifndef COOLFluiD_Physics_NavierStokes_Euler3DRhovt_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DRhovt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/Euler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 3D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class Euler3DRhovt : public Euler3DVarSet {
public: // classes

  /**
   * Constructor
   * @see Euler3D
   */
  Euler3DRhovt(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Euler3DRhovt();
  
  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;
  
  /**
   * Set the jacobians
   */
  virtual void computeJacobians();
  
  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
			     RealMatrix& jacobMin,
			     RealVector& eValues,
			     const RealVector& normal);

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
					 RealMatrix& leftEv,
					 RealVector& eValues,
					 const RealVector& normal);
  
  /**
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const;
  
  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
                            RealVector& result);
  
  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
				     RealVector& result);
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  virtual void computePhysicalData(const Framework::State& state, RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data, Framework::State& state);
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar)
  {
    throw Common::NotImplementedException (FromHere(),"Euler3DRhovt::computePerturbedPhysicalData()");
  }
   
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3; 
  }
  
}; // end of class Euler3DRhovt

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DRhovt_hh
