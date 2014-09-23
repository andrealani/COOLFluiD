#ifndef COOLFluiD_Physics_NEQ_Euler3DNEQPvtyToCons_hh
#define COOLFluiD_Physics_NEQ_Euler3DNEQPvtyToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////
      
/**
 * This class represents a transformer of variables from primitive
 * [p v T] to conservative variables with chemical non equilibrium
 *
 * @author Andrea Lani
 * @author Janos Molnar
 *
 */
class Euler3DNEQPvtyToCons : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerNEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler3DNEQPvtyToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DNEQPvtyToCons();
  
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

  /// acquaintance of the model
  Common::SafePtr<EulerNEQTerm> _model;

  /// array to store density, enthalpy and energy
  RealVector _dhe;

  /// Vector storing the elemental composition
  RealVector _ye;

}; // end of class Euler3DNEQPvtyToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler3DNEQPvtyToCons_hh
