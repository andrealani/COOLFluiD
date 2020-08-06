#ifndef COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaConsToPuvtInPuvt_hh
#define COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaConsToPuvtInPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative [rho rhoU rhoV rhoE rhoK rhoOmega rhoGamma rhoAlpha]
 * to primitive [p u v T K Omega Gamma alpha] starting from primitive variables
 *
 * @author Ray Vandenhoeck
 *
 */
class Euler2DGammaAlphaConsToPuvtInPuvt : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGammaAlphaTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DGammaAlphaConsToPuvtInPuvt
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DGammaAlphaConsToPuvtInPuvt();

  /**
   * Set the transformation matrix from a given state
   */
  void setMatrix(const RealVector& state);

private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

private: //data

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerGammaAlphaTerm> _model;

}; // end of class Euler2DGammaAlphaConsToPuvtInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaConsToPuvtInPuvt_hh
