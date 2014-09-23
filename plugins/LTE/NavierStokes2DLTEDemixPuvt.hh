#ifndef COOLFluiD_Physics_LTE_NavierStokes2DLTEDemixPuvt_hh
#define COOLFluiD_Physics_LTE_NavierStokes2DLTEDemixPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesLTEDemixVarSet.hh"
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
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and LTE Demixing
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
template <class BASE>
class NavierStokes2DLTEDemixPuvt : public NavierStokesLTEDemixVarSet {
public: // classes

  typedef Framework::MultiScalarTerm<BASE> EulerLTEDemixTerm;

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DLTEDemixPuvt(const std::string& name,
                             Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes2DLTEDemixPuvt();

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
   * @pre it is assumed that the composition has already been set
   */
  CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients);

  /**
   * Get the adimensional density
   * @pre the composition will be set here
   */
  CFreal getDensity(const RealVector& state);

  /**
   * Get the adimensional thermal conductivity
   */
  CFreal getThermConductivity(const RealVector& state,
                              const CFreal& dynViscosity);

  /**
   * Get the current update diffusive coefficient
   * @param iEqSS equation subsystem index
   */
  virtual CFreal getCurrUpdateDiffCoeff(CFuint iEqSS);

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
  Common::SafePtr<EulerLTEDemixTerm > _eulerModel;
  
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

  /// array to store the elemental heat transfer coefficients (back up)
  RealVector _lambdaELBkp;

  /// matrix to store the elemental multicomponent difusion coefficients
  /// times mixture density (back up)
  RealMatrix _eldifcoefBkp;

  /// array to store the elemental thermal demixing coefficients
  /// times mixture density (back up)
  RealVector _eltdifcoefBkp;

}; // end of class NavierStokes2DLTEDemixPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DLTEDemixPuvt.ci"

//////////////////////////////////////////////////////////////////////////////
#endif // COOLFluiD_Physics_LTE_NavierStokes2DLTEDemixPuvt_hh
