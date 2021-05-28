#ifndef COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOBSLPuvt_hh
#define COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOBSLPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes3DGReKLogOPuvt.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for primitive
   * variables and K-LogOmega-Gamma-Retheta BSL turbulence model
   *
   * @author Khalil Bensassi 
   * @author Ray Vandenhoeck
   */

class NavierStokes3DGReKLogOBSLPuvt : public NavierStokes3DGReKLogOPuvt {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGReKLogOTerm;

  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DGReKLogOBSLPuvt(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DGReKLogOBSLPuvt();

  /**
   * Compute the blending coeficients needed for BSL and SST variants of the model
   */
  virtual void computeBlendingCoefFromGradientVars(const RealVector& state, const RealVector& gradK, const RealVector& gradOmega);

}; // end of class NavierStokes3DGReKLogOBSLPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_NavierStokes3DGReKLogOBSLPuvt_hh
