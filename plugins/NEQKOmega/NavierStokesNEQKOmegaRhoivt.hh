#ifndef COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmegaRhoivt_hh
#define COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmegaRhoivt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NEQKOmega/NavierStokesNEQKOmegaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

      namespace NEQKOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model for thermo and/or chemical
   * NEQ with K-Omega turbulence model
   *
   * @author Andrea Lani
   */
template <typename BASE>      
class NavierStokesNEQKOmegaRhoivt : public NavierStokesNEQKOmegaVarSet<BASE> {
public: // classes
  
  /**
   * Constructor
   */
  NavierStokesNEQKOmegaRhoivt(const std::string& name,
			      Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesNEQKOmegaRhoivt();
  
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

  /**
   * Set velocity, pressure, temperature, density variables
   * starting from state variables
   */
  virtual void setGradientState(const RealVector& state);
  
protected:

  /// vibrational temperatures
  RealVector _tempVib;
    
}; // end of class NavierStokesNEQKOmegaRhoivt
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace   namespace NEQKOmega {

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQKOmegaRhoivt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQKOmega_NavierStokesNEQKOmegaRhoivt_hh
