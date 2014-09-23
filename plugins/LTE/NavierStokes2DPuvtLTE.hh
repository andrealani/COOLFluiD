#ifndef COOLFluiD_Physics_LTE_NavierStokes2DPuvtLTE_hh
#define COOLFluiD_Physics_LTE_NavierStokes2DPuvtLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
template <typename CTERM>
class NavierStokes2DPuvtLTE : public NavierStokes::NavierStokes2DVarSet {
public: // classes

  typedef CTERM CONVTERM;

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DPuvtLTE(const std::string& name,
			Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokes2DPuvtLTE();

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

  /**
   * Get the adimensional dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients);

  /**
   * Get the adimensional density
   * @pre the composition will be set here
   */
  virtual CFreal getDensity(const RealVector& state);

  /**
   * Get the adimensional thermal conductivity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getThermConductivity(const RealVector& state,
			      const CFreal& dynViscosity);

  /**
   * set up private data
   */
  virtual void setup();

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

private :

  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// convective model
  Common::SafePtr<CTERM> _eulerModel;

  /// array for the composition
  RealVector _tempX;

}; // end of class NavierStokes2DPuvtLTE

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DPuvtLTE.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_NavierStokes2DPuvtLTE_hh
