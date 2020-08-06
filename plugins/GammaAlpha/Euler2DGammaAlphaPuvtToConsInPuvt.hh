#ifndef COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaPuvtToConsInPuvt_hh
#define COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaPuvtToConsInPuvt_hh

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
 * This class represents a transformer of variables from primitive [p u v T K Omega Gamma alpha]
 * to conservative [rho rhoU rhoV rhoE rhoK rhoOmega rhoGamma rhoAlpha]  starting from primitive variables
 *
 * @author Ray Vandenhoeck
 *
 */
class Euler2DGammaAlphaPuvtToConsInPuvt : public Framework::VarSetMatrixTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerGammaAlphaTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DGammaAlphaPuvtToConsInPuvt
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DGammaAlphaPuvtToConsInPuvt();

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

}; // end of class Euler2DGammaAlphaPuvtToConsInPuvt

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_Euler2DGammaAlphaPuvtToConsInPuvt_hh
