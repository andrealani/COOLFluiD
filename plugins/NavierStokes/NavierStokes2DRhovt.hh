#ifndef COOLFluiD_Physics_NavierStokes_NavierStokes2DRhovt_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokes2DRhovt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {
      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class NavierStokes2DRhovt : public NavierStokes2DVarSet {
public:

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DRhovt(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes2DRhovt();

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  void setComposition(const RealVector& state,
		      const bool isPerturb,
		      const CFuint iVar)
  {
  }

  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  void setGradientVars(const std::vector<RealVector*>& states,
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

  /**
   * Get the adimensional dynamic viscosity
   */
  CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients);

  /**
   * Get the adimensional density
   */
  CFreal getDensity(const RealVector& state);

  /**
   * Get the adimensional thermal conductivity
   */
  CFreal getThermConductivity(const RealVector& state,
			      const CFreal& dynViscosity)
  {
    if (Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional()) {
      return dynViscosity;
    }
    return dynViscosity*getModel().getCpOverPrandtl();
  }

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

private:

  /// convective model
  Common::SafePtr<EulerTerm> _eulerModel;

}; // end of class NavierStokes2DRhovt

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokes2DRhovt_hh
