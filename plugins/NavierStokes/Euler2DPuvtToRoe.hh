#ifndef COOLFluiD_Physics_NavierStokes_Euler2DPuvtToRoe_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DPuvtToRoe_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from {p u v T}
 * to Roe variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DPuvtToRoe : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DPuvtToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DPuvtToRoe();
  
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
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler2DPuvtToRoe

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DPuvtToRoe_hh
