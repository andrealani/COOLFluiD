#ifndef COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaSSTPuvt_hh
#define COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaSSTPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes3DGammaAlphaBSLPuvt.hh"
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

class NavierStokes3DGammaAlphaSSTPuvt : public NavierStokes3DGammaAlphaBSLPuvt {
public: // classes

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGammaAlphaTerm;

  /**
   * Constructor
   * @see NavierStokes3D
   */
  NavierStokes3DGammaAlphaSSTPuvt(const std::string& name,
			     Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~NavierStokes3DGammaAlphaSSTPuvt();

  /**
   * Set the coefficients for the model
   */
  virtual void setModelCoefficients();

  /**
   * Get the adimensional turbulent dynamic viscosity
   * @pre it is assumed that the composition has already been set
   */
  virtual CFreal getTurbDynViscosityFromGradientVars(const RealVector& state, const std::vector<RealVector*>& gradients);

}; // end of class NavierStokes3DGammaAlphaSSTPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_NavierStokes3DGammaAlphaSSTPuvt_hh
