#ifndef COOLFluiD_Physics_GReKO_Euler2DGReKOConsToRoe_hh
#define COOLFluiD_Physics_GReKO_Euler2DGReKOConsToRoe_hh

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
 * This class represents a transformer of variables from conservative
 * to Roe variables for Gamma-ReTheta Transition model
 *
 * @author Khalil Bensassi
 *
 */

class Euler2DGReKOConsToRoe : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerGReKOTerm;

  /**
   * Default constructor without arguments
   */
  Euler2DGReKOConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DGReKOConsToRoe();
  
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
  Common::SafePtr<EulerGReKOTerm> _model;

}; // end of class Euler2DGReKOConsToRoe

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GReKO_Euler2DGReKOConsToRoe_hh
