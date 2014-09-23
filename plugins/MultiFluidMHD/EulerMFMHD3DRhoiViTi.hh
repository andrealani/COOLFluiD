#ifndef COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD3DRhoiViTi_hh
#define COOLFluiD_Physics_MultiFluidMHD_EulerMFMHD3DRhoiViTi_hh

//////////////////////////////////////////////////////////////////////////////

#include "MultiFluidMHDVarSet.hh"
#include "Maxwell/Maxwell3DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a multi-fluid Euler physical 
 * model 3D for Rhoi, Vi, Ti
 * variables:
 * Bx, By, Bz, Ex, Ey, Ez, Psi, Phi,
 * {rhoi}, {Ui, Vi, Wi}, {Ti}
 * 
 * PHYSICAL DATA ARRAY:
 * Bx, By, Bz, Ex, Ey, Ez, Psi, Phi,
 * rho, XP, YP, ZP,
 * {yi},
 * {Ui, Vi, Wi},
 * {Ti, pi, ai, Hi} 
 * 
 * 
 * @author Alejandro Alvarez
 * 
 */
class EulerMFMHD3DRhoiViTi : public MultiFluidMHDVarSet<Maxwell::Maxwell3DProjectionVarSet> {
public: // classes

  typedef Framework::MultiScalarTerm<EulerMFMHDTerm> PTERM;
  
  /**
   * Constructor
   * @see MultiFluidMHDModel
   */
  EulerMFMHD3DRhoiViTi(Common::SafePtr<Framework::BaseTerm> term);
  
  /**
   * Default destructor
   */
  virtual ~EulerMFMHD3DRhoiViTi();
  
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
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const;
  
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
   * @see EulerPhysicalModel
   */
  void computePhysicalData(const Framework::State& state, RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
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
  
  /// Array with the particle mass of the species
  RealVector			   _m_i;
  
}; // end of class EulerMFMHD3DRhoiViTi
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHD3DRhoiViTi_hh
