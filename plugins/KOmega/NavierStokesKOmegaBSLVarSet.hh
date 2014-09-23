#ifndef COOLFluiD_Physics_KOmega_NavierStokesKOmegaBSLVarSet_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKOmegaBSLVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/NavierStokesKOmegaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and K-Omega BSL turbulence model
   *
   * @author Thomas Wuilbaut
   */
template <typename BASE, int SGROUP>
class NavierStokesKOmegaBSLVarSet : public NavierStokesKOmegaVarSet<BASE, SGROUP> {
public: // classes

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokesKOmegaBSLVarSet(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokesKOmegaBSLVarSet();
  
  /**
   * Compute the blending coeficients needed for BSL and SST variants of the model
   */
  virtual void computeBlendingCoefFromGradientVars(const RealVector& state, 
						   const RealVector& gradK, 
						   const RealVector& gradOmega);
  
}; // end of class NavierStokesKOmegaBSLVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKOmegaBSLVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKOmegaBSLVarSet_hh
