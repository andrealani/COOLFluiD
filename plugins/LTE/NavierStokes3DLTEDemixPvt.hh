#ifndef COOLFluiD_Physics_LTE_NavierStokes3DLTEDemixPvt_hh
#define COOLFluiD_Physics_LTE_NavierStokes3DLTEDemixPvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesLTEDemixVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables and LTE Demixing
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
class NavierStokes3DLTEDemixPvt : public NavierStokesLTEDemixVarSet {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerLTEDemixTerm;

  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DLTEDemixPvt(const std::string& name,
                            Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DLTEDemixPvt();

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
   * set up private data
   */
  void setup();

  /**
   * Get the diffusive flux
   */
  RealVector& getFlux(const RealVector& values,
                      const std::vector<RealVector*>& gradients,
                      const RealVector& normal,
                      const CFreal& radius);

  /**
   * Get the diffusive flux vector
   */
  RealMatrix& getFlux(const RealVector& values,
                      const std::vector<RealVector*>& gradients,
                      const CFreal& radius);

  /**
   * Get the heat flux
   */
  CFreal getHeatFlux(const RealVector& state,
		     const std::vector<RealVector*>& gradients,
		     const RealVector& normal);

protected:

  /**
   * Set the gradient variables starting from state variables
   */
  virtual void setGradientState(const RealVector& state);

private :

  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  /// convective model
  Common::SafePtr<EulerLTEDemixTerm> _eulerModel;

  /// array for the composition
  RealVector _tempX;

  /// array to store the mass fractions of elements
  RealVector _ye;

  /// array to store the elemental heat transfer coefficients
  RealVector _lambdaEL;

  /// matrix to store the elemental multicomponent difussion coefficients
  /// times mixture density
  RealMatrix _eldifcoef;

  /// array to store the elemental thermal demixing coefficients
  /// times mixture density
  RealVector _eltdifcoef;

}; // end of class NavierStokes3DLTEDemixPvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_NavierStokes3DLTEDemixPvt_hh
