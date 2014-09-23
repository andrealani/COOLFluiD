#ifndef COOLFluiD_Physics_NavierStokes_Euler3DConsToPrim_hh
#define COOLFluiD_Physics_NavierStokes_Euler3DConsToPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

      class EulerTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to primitive [rho u v w p] variables
 *
 * @author Willem Deconinck
 *
 */
class Euler3DConsToPrim : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler3DConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DConsToPrim();

  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);

  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& pdata, Framework::State& result);
  
private:

  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerTerm> _model;

}; // end of class Euler3DConsToPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler3DConsToPrim_hh
