#ifndef COOLFluiD_Physics_GReKO_Euler3DGReKLogOPuvtToCons_hh
#define COOLFluiD_Physics_GReKO_Euler3DGReKLogOPuvtToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [p u v w T K Omega] to conservative variables [rho rhoU rhoV rhoE rhoK rhoOmega]
 *
 * @author Khalil Bensassi
 * @author Ray Vandenhoeck
 *
 */
class Euler3DGReKLogOPuvtToCons : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerGReKLogOTerm;

  /**
   * Constructor
   */
  Euler3DGReKLogOPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DGReKLogOPuvtToCons();
  
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
  Common::SafePtr<EulerGReKLogOTerm> _model;

}; // end of class Euler3DGReKLogOPuvtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_Euler3DGReKLogOPuvtToCons_hh
