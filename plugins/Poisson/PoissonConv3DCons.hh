#ifndef COOLFluiD_Physics_Poisson_PoissonConv3DCons_hh
#define COOLFluiD_Physics_Poisson_PoissonConv3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "PoissonConv3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Poisson physical model 3D for conservative
 * variables
 *
 * @author Alejandro Alvarez Laguna
 */
class PoissonConv3DCons : public PoissonConv3DVarSet {
public: // classes

  /**
   * Constructor
   * @see PoissonModel
   */
  PoissonConv3DCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~PoissonConv3DCons();

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
   * Give dimensional values to the adimensional state variables
   */
  void setDimensionalValues(const Framework::State& state, RealVector& result);
  
  /**
   * Give adimensional values to the dimensional state variables
   */
  void setAdimensionalValues(const Framework::State& state, RealVector& result);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see PoissonPhysicalModel
   */
  void computePhysicalData(const Framework::State& state, RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see PoissonPhysicalModel
   */
  void computeStateFromPhysicalData(const RealVector& data, Framework::State& state);
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar);
  
  /**
   * Returns true if the state fed doesn't have unphysical values
   * Overrides the standard definition in ConvectiveVarSet.
   */
  bool isValid(const RealVector& data);
  
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
  
}; // end of class PoissonConv3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonConv3DCons_hh
