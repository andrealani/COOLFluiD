#ifndef COOLFluiD_Physics_KOmega_Euler3DKOmegaConsToRoe_hh
#define COOLFluiD_Physics_KOmega_Euler3DKOmegaConsToRoe_hh

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
 * This class represents a transformer of variables from conservative
 * to Roe variables
 *
 * @author Thomas Wuilbaut
 * @author Milan Zaloudek
 *
 */

class Euler3DKOmegaConsToRoe : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> EulerKOmegaTerm;

  /**
   * Default constructor without arguments
   */
  Euler3DKOmegaConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DKOmegaConsToRoe();
  
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
  Common::SafePtr<EulerKOmegaTerm> _model;

}; // end of class Euler3DKOmegaConsToRoe

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_Euler3DKOmegaConsToRoe_hh
