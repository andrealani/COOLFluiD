#ifndef COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaPvt_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model  for primitive
   * variables and K-LogOmega turbulence model
   *
   * @author Thomas Wuilbaut
   * @author Khalil Bensassi
   * @author Ray Vandenhoeck
   */
template <typename BASE>      
class NavierStokesKLogOmegaPvt : public BASE {
public: // classes
  
  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesKLogOmegaPvt(const std::string& name,
			Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesKLogOmegaPvt();
  
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
  
  /**
   * Get the adimensional laminar dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getLaminarDynViscosityFromGradientVars(const RealVector& state);
  
  
  /**
   * Get the adimensional density
   * @pre the composition will be set here
   */
  virtual CFreal getDensity(const RealVector& state);
  
  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
				      const CFreal& dynViscosity)
  {
    if (Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional()) {
      return dynViscosity;
    }
    return dynViscosity*this->getModel().getCpOverPrandtl();
  }
  
  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);
  
  /**
   * set up private data
   */
  virtual void setup();
  
  /// Compute the speed of sound
  virtual CFreal getSqSoundSpeed(const RealVector& state);
  
protected:
  
  /// Set the gradient variables starting from state variables
  virtual void setGradientState(const RealVector& state);
  
}; // end of class NavierStokesKLogOmegaPvt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKLogOmegaPvt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaPvt_hh
