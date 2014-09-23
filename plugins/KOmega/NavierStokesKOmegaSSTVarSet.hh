#ifndef COOLFluiD_Physics_KOmega_NavierStokesKOmegaSSTVarSet_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKOmegaSSTVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/NavierStokesKOmegaBSLVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model with Menter's 
   * K-Omega SST turbulence model
   *
   * @author Thomas Wuilbaut
   */
template <typename BASE, int SGROUP>
class NavierStokesKOmegaSSTVarSet : public NavierStokesKOmegaBSLVarSet<BASE, SGROUP> {
public: // classes

  /**
   * Constructor
   */
  NavierStokesKOmegaSSTVarSet(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokesKOmegaSSTVarSet();
  
  /**
   * Set the coefficients for the model
   */
  virtual void setModelCoefficients();

  /**
   * Get the adimensional turbulent dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getTurbDynViscosityFromGradientVars(const RealVector& state, 
						     const std::vector<RealVector*>& gradients);
private:
  
  /// square
  CFreal sq(CFreal a) const {return a*a;}
  
}; // end of class NavierStokesKOmegaSSTVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKOmegaSSTVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKOmegaSSTVarSet_hh
