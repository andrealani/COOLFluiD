#ifndef COOLFluiD_Physics_LTE_NavierStokes3DPvtLTE_hh
#define COOLFluiD_Physics_LTE_NavierStokes3DPvtLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokes3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Physics {

    namespace NavierStokes {
      class EulerTerm;
    }

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class NavierStokes3DPvtLTE : public NavierStokes::NavierStokes3DVarSet {
public: // classes

  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DPvtLTE(const std::string& name,
		       Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DPvtLTE();

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  void setComposition(const RealVector& state,
		      const bool isPerturb,
		      const CFuint iVar);

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
			      const CFreal& dynViscosity);

  /**
   * set up private data
   */
  void setup();

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

private :

  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// convective model
  Common::SafePtr<NavierStokes::EulerTerm> _eulerModel;

  /// array for the composition
  RealVector _tempX;

}; // end of class NavierStokes3DPvtLTE

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_NavierStokes3DPvtLTE_hh
