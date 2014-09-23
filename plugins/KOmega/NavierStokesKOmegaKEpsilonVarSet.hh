#ifndef COOLFluiD_Physics_KOmega_NavierStokesKOmegaKEpsilonVarSet_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKOmegaKEpsilonVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/NavierStokesKOmegaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and K-Omega KEpsilon turbulence model
   *
   * @author Thomas Wuilbaut
   */
template <typename BASE, int SGROUP>
class NavierStokesKOmegaKEpsilonVarSet : public NavierStokesKOmegaVarSet<BASE, SGROUP> {
public: // classes

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokesKOmegaKEpsilonVarSet(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokesKOmegaKEpsilonVarSet();
  
  /**
   * Compute the blending coeficients needed for BSL and SST variants of the model
   */
  virtual void computeBlendingCoefFromGradientVars(const RealVector& state, 
						   const RealVector& gradK, 
						   const RealVector& gradOmega, 
						   const CFreal& distance);
  
}; // end of class NavierStokesKOmegaKEpsilonVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKOmegaKEpsilonVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKOmegaKEpsilonVarSet_hh
