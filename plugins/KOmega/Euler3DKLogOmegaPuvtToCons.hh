#ifndef COOLFluiD_Physics_KOmega_Euler3DKLogOmegaPuvtToCons_hh
#define COOLFluiD_Physics_KOmega_Euler3DKLogOmegaPuvtToCons_hh

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
 * [p u v w T K Omega] to conservative variables [rho rhoU rhoV rhoW rhoE rhoK rhoOmega]
 *
 * @author Thomas Wuilbaut
 * @author Milan Zaloudek
 * @author Ray Vandenhoeck
 *
 */
class Euler3DKLogOmegaPuvtToCons : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerKLogOmegaTerm;

  /**
   * Constructor
   */
  Euler3DKLogOmegaPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DKLogOmegaPuvtToCons();
  
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

}; // end of class Euler3DKLogOmegaPuvtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace KLogOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_Euler3DKLogOmegaPuvtToCons_hh
