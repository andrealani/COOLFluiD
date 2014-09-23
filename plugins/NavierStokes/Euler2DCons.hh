#ifndef COOLFluiD_Physics_NavierStokes_Euler2DCons_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative
 * variables
 *
 * @author Andrea Lani
 */
class Euler2DCons : public Euler2DVarSet {
public: // classes

  /**
   * Constructor
   * @see Euler2D
   */
  Euler2DCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~Euler2DCons();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const;
  
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
   * Set the first right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  void setEigenVect1(RealVector& r1, Framework::State& state,const RealVector& normal);
  
  /**
   * Set the second right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}\f$)
   */
  void setEigenVect2(RealVector& r2, Framework::State& state, const RealVector& normal);
  
  /**
   * Set the third right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  void setEigenVect3(RealVector& r3, Framework::State& state, const RealVector& normal);
  
  /**
   * Set the fourth right eigen vector (corresponding to \f$\vec{u} \cdot \vec{n}+ a \f$)
   */
  void setEigenVect4(RealVector& r4, Framework::State& state, const RealVector& normal);
  
  /**
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const;
  
  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state, RealVector& result);
  
  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state, RealVector& result);
  
  
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
					    CFuint iVar);
  
  /**
   * Returns true if the state fed doesn't have unphysical values
   * Overrides the standard definition in ConvectiveVarSet.
   */
  virtual bool isValid(const RealVector& data);

  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(2); velIDs[XX] = 1; velIDs[YY] = 2;
  }
  
protected:
  
  /**
   * Set the constant part (independent form the solution) of the
   * jacobians
   */
  virtual void setConstJacob();
  
private: // data

  /// temporary matrix of right eigenvalues
  RealMatrix                       _rightEv;

  /// temporary matrix of left eigenvalues
  RealMatrix                       _leftEv;

}; // end of class Euler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DCons_hh
