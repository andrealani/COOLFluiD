#ifndef COOLFluiD_Physics_KOmega_Euler2DKLogOmegaPuvtToCons_hh
#define COOLFluiD_Physics_KOmega_Euler2DKLogOmegaPuvtToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [p u v T K LogOmega] to conservative variables [rho rhoU rhoV rhoE rhoK rhoLogOmega]
 *
 * @author Thomas Wuilbaut
 * @author Ray Vandenhoeck
 *
 */
class Euler2DKLogOmegaPuvtToCons : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerKLogOmegaTerm;

  /**
   * Constructor
   */
  Euler2DKLogOmegaPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DKLogOmegaPuvtToCons();
  
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
  Common::SafePtr<EulerKLogOmegaTerm> _model;

}; // end of class Euler2DKLogOmegaPuvtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_Euler2DKLogOmegaPuvtToCons_hh
