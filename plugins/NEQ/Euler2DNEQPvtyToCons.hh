#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQPvtyToCons_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQPvtyToCons_hh

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
 * [p u v T] to conservative variables with LTE and Demixing
 *
 * @author Andrea Lani
 * @author Janos Molnar
 *
 */
class Euler2DNEQPvtyToCons : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerNEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQPvtyToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQPvtyToCons();
  
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

}; // end of class Euler2DNEQPvtyToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQPvtyToCons_hh
