#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtToCons_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtToCons_hh

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
 *
 */
class Euler2DNEQRhoivtToCons : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerNEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQRhoivtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQRhoivtToCons();

  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
protected:

  /// acquaintance of the model
  Common::SafePtr<EulerNEQTerm> _model;
 
  /// array to store density, enthalpy and energy
  RealVector _dhe;

  /// Vector storing the elemental composition
  RealVector _ye;

}; // end of class Euler2DNEQRhoivtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtToCons_hh
