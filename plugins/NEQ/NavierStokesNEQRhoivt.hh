#ifndef COOLFluiD_Physics_NEQ_NavierStokesNEQRhoivt_hh
#define COOLFluiD_Physics_NEQ_NavierStokesNEQRhoivt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model  for primitive
   * variables and chemical NEQ
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
template <class BASE>
class NavierStokesNEQRhoivt : public BASE {
public: // classes

  /**
   * Constructor
   */
  NavierStokesNEQRhoivt(const std::string& name,
			  Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokesNEQRhoivt();

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
  void setGradientVarGradients(const std::vector<RealVector*>& states,
                               const std::vector< std::vector<RealVector*> >& stateGradients,
                               std::vector< std::vector<RealVector*> >& gradVarGradients,
                               const CFuint stateSize);

  /**
   * Compute the gradients of the states starting from gradient variable gradients (pressure, velocity, temperature)
   */
  void setStateGradients(const std::vector<RealVector*>& states,
                         const std::vector< std::vector<RealVector*> >& gradVarGradients,
                         std::vector< std::vector<RealVector*> >& stateGradients,
                         const CFuint stateSize);
  
  /// Get the diffusive flux jacobian
  void computeFluxJacobian(const RealVector& state,
			   const RealVector& gradientJacob,
			   const RealVector& normal,
			   const CFreal& radius,
			   RealMatrix& fluxJacob);
  
  /**
   * Get the adimensional dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getDynViscosity(const RealVector& state,
				 const std::vector<RealVector*>& gradients);

  /**
   * Get the adimensional density
   * @pre the composition will be set here
   */
  virtual CFreal getDensity(const RealVector& state);

  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
				      const CFreal& dynViscosity);

  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);

  /**
   * set up private data
   */
  virtual void setup();

protected:

  /**
   * Set velocity, pressure, temperature, density variables
   * starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

protected:

  /// vibrational temperatures
  RealVector _tempVib;

}; // end of class NavierStokesNEQRhoivt

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQRhoivt.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_NavierStokesNEQRhoivt_hh
