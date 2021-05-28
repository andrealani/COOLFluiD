#ifndef COOLFluiD_Physics_GReKO_Euler3DGReKLogOPuvt_hh
#define COOLFluiD_Physics_GReKO_Euler3DGReKLogOPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/Euler3DPvt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 3D for primitive
   * variables with Gamma-Re transition Model coupled with k-LogOmega turbulence model
   *
   * @author Khalil Bensassi
   * @author Ray Vandenhoeck
   */

class Euler3DGReKLogOPuvt : public Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler3DPvt<Physics::NavierStokes::Euler3DVarSet> > {
public: // classes

  /**
   * Puvttructor
   * @see Euler3D
   */
  Euler3DGReKLogOPuvt(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Euler3DGReKLogOPuvt();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  void setup();

  /**
   * Get extra variable names
   */
  std::vector<std::string> getExtraVarNames() const;

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;

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
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result);

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result);


  /**
   * Set other adimensional values for useful physical quantities
   */
  void setDimensionalValuesPlusExtraValues
  (const Framework::State& state, RealVector& result,
   RealVector& extra);
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar);
  
  /**
   * Set the PhysicalData corresponding to the given State
   * @see EulerPhysicalModel
   */
  virtual void computePhysicalData(const Framework::State& state,
			       RealVector& data);
  
  /**
   * Set a State starting from the given PhysicalData
   * @see EulerPhysicalModel
   */
  virtual void computeStateFromPhysicalData(const RealVector& data,
				       Framework::State& state);
     
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  {
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3;
  }
  
}; // end of class Euler3DGReKLogOPuvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_Euler3DGReKLogOPuvt_hh
