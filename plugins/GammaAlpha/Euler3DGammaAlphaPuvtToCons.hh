#ifndef COOLFluiD_Physics_GammaAlpha_Euler3DGammaAlphaPuvtToCons_hh
#define COOLFluiD_Physics_GammaAlpha_Euler3DGammaAlphaPuvtToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [p u v T K Omega gamma alpha] to conservative variables [rho rhoU rhoV rhoE rhoK rhoOmega rhoGamma rhoAlpha]
 *
 * @author Ray Vandenhoeck
 *
 */
class Euler3DGammaAlphaPuvtToCons : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerGammaAlphaTerm;

  /**
   * Constructor
   */
  Euler3DGammaAlphaPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DGammaAlphaPuvtToCons();
  
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerGammaAlphaTerm> _model;

}; // end of class Euler3DGammaAlphaPuvtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_Euler3DGammaAlphaPuvtToCons_hh
