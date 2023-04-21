#ifndef COOLFluiD_Physics_MHD_MHD3DProjectionDiffPrim_hh
#define COOLFluiD_Physics_MHD_MHD3DProjectionDiffPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD3DProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {
      class MHDProjectionTerm;
      
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a MHD physical model 3D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class MHD3DProjectionDiffPrim : public MHD3DProjectionDiffVarSet {
public: // classes

  /**
   * Constructor
   * @see MHD3DProjectionDiffVarSet
   */
  MHD3DProjectionDiffPrim(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~MHD3DProjectionDiffPrim();

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  virtual void setComposition(const RealVector& state,
			      const bool isPerturb,
			      const CFuint iVar);
  
  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  virtual void setGradientVars(const std::vector<RealVector*>& states,
			       RealMatrix& values,
			       const CFuint stateSize);
  
  /**
   * Compute required gradients (velocity, Temperature) starting from the gradients of the states
   */
  virtual void setGradientVarGradients(const std::vector<RealVector*>& states,
				       const std::vector< std::vector<RealVector*> >& stateGradients,
				       std::vector< std::vector<RealVector*> >& gradVarGradients,
				       const CFuint stateSize);

  /**
   * Compute the gradients of the states starting from gradient variable gradients (pressure, velocity, temperature)
   */
  virtual void setStateGradients(const std::vector<RealVector*>& states,
				 const std::vector< std::vector<RealVector*> >& gradVarGradients,
				 std::vector< std::vector<RealVector*> >& stateGradients,
				 const CFuint stateSize);

  /// Get the diffusive flux
  virtual void computeFluxJacobian(const RealVector& state,
				   const RealVector& gradientJacob,
				   const RealVector& normal,
				   const CFreal& radius,
				   RealMatrix& fluxJacob);
  
  /**
   * Get the adimensional dynamic viscosity
   */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients);

  /**
   * Get the adimensional density
   */
  virtual CFreal getDensity(const RealVector& state);

  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
				      const CFreal& dynViscosity)
  {
  }

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

private:

  /// convective model
  Common::SafePtr<MHDProjectionTerm> _eulerModel;
  
  /// array for composition
  RealVector _tempX;
  
}; // end of class MHD3DProjectionDiffPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHD3DProjectionDiffPrim_hh
