#ifndef COOLFluiD_Physics_LTE_NavierStokesLTEDemixVarSet_hh
#define COOLFluiD_Physics_LTE_NavierStokesLTEDemixVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model for primitive
   * variables and LTE Demixing
   *
   * @author Andrea Lani
   * @author Janos Molnar
   */
class NavierStokesLTEDemixVarSet : public NavierStokes::NavierStokesVarSet {
public: // classes

  /**
   * Constructor
   * @see NavierStokes
   */
  NavierStokesLTEDemixVarSet(const std::string& name,
	  Common::SafePtr<Framework::PhysicalModelImpl> model) :
    NavierStokes::NavierStokesVarSet(name, model)
  {
  }

  /**
   * Default destructor
   */
  virtual ~NavierStokesLTEDemixVarSet()
  {
  }

  /**
   * Set up private data
   */
  virtual void setup()
  {
  }

  /**
   * Set the composition
   * @pre this function has to be called before any other function
   *      computing other physical quantities
   */
  virtual void setComposition(const RealVector& state,
			      const bool isPerturb,
			      const CFuint iVar) = 0;

  /**
   * Set the quantities needed to compute gradients (pressure,
   * velocity, etc.) starting from the states
   */
  virtual void setGradientVars(const std::vector<RealVector*>& states,
			       RealMatrix& values,
			       const CFuint stateSize) = 0;

  /**
   * Get the adimensional dynamic viscosity
   */
  virtual CFreal getDynViscosity(const RealVector& state, const std::vector<RealVector*>& gradients) = 0;

  /**
   * Get the adimensional density
   */
  virtual CFreal getDensity(const RealVector& state) = 0;

  /**
   * Get the adimensional thermal conductivity
   */
  virtual CFreal getThermConductivity(const RealVector& state,
				      const CFreal& dynViscosity)
  {
    throw Common::NotImplementedException
      (FromHere(), "NavierStokesLTEDemixVarSet::getThermConductivity()");
  }

  /**
   * Get the diffusive flux
   */
  virtual RealVector& getFlux(const RealVector& values,
			      const std::vector<RealVector*>& gradients,
			      const RealVector& normal,
			      const CFreal& radius) = 0;

}; // end of class NavierStokesLTEDemixVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_NavierStokesLTEDemixVarSet_hh
