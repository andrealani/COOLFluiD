#ifndef COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaBSLPuvt_hh
#define COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaBSLPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes3DGammaAlphaPuvt.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NavierStokes physical model 3D for gamma-alpha
   *
   * @author Ray Vandenhoeck
   */

class NavierStokes3DGammaAlphaBSLPuvt : public NavierStokes3DGammaAlphaPuvt {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGammaAlphaTerm;

  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DGammaAlphaBSLPuvt(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DGammaAlphaBSLPuvt();

  /**
   * Compute the blending coeficients needed for BSL and SST variants of the model
   */
  virtual void computeBlendingCoefFromGradientVars(const RealVector& state, const RealVector& gradK, const RealVector& gradOmega);
  
  /**
   * Get the adimensional turbulent dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getTurbDynViscosityFromGradientVars(const RealVector& state, const std::vector<RealVector*>& gradients);

}; // end of class NavierStokes3DGammaAlphaBSLPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaBSLPuvt_hh
