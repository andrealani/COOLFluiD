#ifndef COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaSSTVarSet_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaSSTVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/NavierStokesKLogOmegaBSLVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model with Menter's 
   * K-LogOmega SST turbulence model
   *
   * @author Thomas Wuilbaut
   * @author Ray Vandenhoeck
   */
template <typename BASE, int SGROUP>
class NavierStokesKLogOmegaSSTVarSet : public NavierStokesKLogOmegaBSLVarSet<BASE, SGROUP> {
public: // classes

  /**
   * Constructor
   */
  NavierStokesKLogOmegaSSTVarSet(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokesKLogOmegaSSTVarSet();
  
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
  
}; // end of class NavierStokesKLogOmegaSSTVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKLogOmegaSSTVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaSSTVarSet_hh
