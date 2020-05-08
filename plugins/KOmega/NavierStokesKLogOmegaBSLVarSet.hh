#ifndef COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaBSLVarSet_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaBSLVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "KOmega/NavierStokesKLogOmegaVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and K-LogOmega BSL turbulence model
   *
   * @author Thomas Wuilbaut
   * @author Ray Vandenhoeck
   */
template <typename BASE, int SGROUP>
class NavierStokesKLogOmegaBSLVarSet : public NavierStokesKLogOmegaVarSet<BASE, SGROUP> {
public: // classes

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokesKLogOmegaBSLVarSet(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~NavierStokesKLogOmegaBSLVarSet();
  
  /**
   * Compute the blending coeficients needed for BSL and SST variants of the model
   */
  virtual void computeBlendingCoefFromGradientVars(const RealVector& state, 
						   const RealVector& gradK, 
						   const RealVector& gradOmega);
  
}; // end of class NavierStokesKLogOmegaBSLVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesKLogOmegaBSLVarSet.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaBSLVarSet_hh
