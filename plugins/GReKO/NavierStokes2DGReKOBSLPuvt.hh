#ifndef COOLFluiD_Physics_GReKO_NavierStokes2DGReKOBSLPuvt_hh
#define COOLFluiD_Physics_GReKO_NavierStokes2DGReKOBSLPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DGReKOPuvt.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 2D for primitive
   * variables and K-Omega-Gamma-Retheta BSL turbulence model
   *
   * @author Khalil Bensassi 
   */

class NavierStokes2DGReKOBSLPuvt : public NavierStokes2DGReKOPuvt {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGReKOTerm;

  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DGReKOBSLPuvt(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes2DGReKOBSLPuvt();

  /**
   * Compute the blending coeficients needed for BSL and SST variants of the model
   */
  virtual void computeBlendingCoefFromGradientVars(const RealVector& state, const RealVector& gradK, const RealVector& gradOmega);

}; // end of class NavierStokes2DGReKOBSLPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokes2DGReKOBSLPuvt_hh
