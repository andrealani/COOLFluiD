#ifndef COOLFluiD_Physics_NavierStokes_Euler2DConsToPrim_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DConsToPrim_hh

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
 * to primitive [rho u v p] variables
 *
 * @author Kris Van den Abeele
 *
 */
class Euler2DConsToPrim : public Framework::VarSetTransformer {
public:

  /**
   * Default constructor without arguments
   */
  Euler2DConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DConsToPrim();

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

}; // end of class Euler2DConsToPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DConsToPrim_hh
