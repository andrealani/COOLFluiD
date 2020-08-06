#ifndef COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaCons_hh
#define COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 2D for conservative
   * variables with Gamma-alpha transition model
   *
   * @author Ray Vandenhoeck
   */
class Euler2DGammaAlphaCons : public Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DCons> {
public: // classes

  /**
   * Constructor
   * @see Euler2D
   */
  Euler2DGammaAlphaCons(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Euler2DGammaAlphaCons();

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
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Set the constant part of the jacobians
   */
  virtual void setConstJacob();

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
    velIDs.resize(2); velIDs[XX] = 1; velIDs[YY] = 2;
  }
  
}; // end of class Euler2DGammaAlphaCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaCons_hh
