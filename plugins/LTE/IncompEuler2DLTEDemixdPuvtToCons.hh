#ifndef COOLFluiD_Physics_LTE_IncompEuler2DLTEDemixdPuvtToCons_hh
#define COOLFluiD_Physics_LTE_IncompEuler2DLTEDemixdPuvtToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/IncompEulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {
      class IncompEulerTerm;
    }

    namespace LTE {

/////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [dp u v T] to conservative variables with LTE and Demixing
 *
 * @author Andrea Lani
 * @author Janos Molnar
 *
 */
class IncompEuler2DLTEDemixdPuvtToCons : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::IncompEulerTerm> IncompEulerLTEDemixTerm;

  /**
   * Default constructor without arguments
   */

  IncompEuler2DLTEDemixdPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~IncompEuler2DLTEDemixdPuvtToCons();

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
  Common::SafePtr<IncompEulerLTEDemixTerm> _model;

  /// array to store density, enthalpy and energy
  RealVector _dhe;

  /// Vector storing the elemental composition
  RealVector _ye;

}; // end of class Euler2DLTEDemixPuvtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_Euler2DLTEDemixPuvtToCons_hh




