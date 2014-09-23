#ifndef COOLFluiD_Physics_SA_Euler2DSAPuvtToCons_hh
#define COOLFluiD_Physics_SA_Euler2DSAPuvtToCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [p u v T K] to conservative variables [rho rhoU rhoV rhoE rhoK]
 *
 * @author Joao Pinto
 * @author Thomas Wuilbaut
 *
 */
class Euler2DSAPuvtToCons : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerSATerm;
  
  /**
   * Constructor
   */
  Euler2DSAPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DSAPuvtToCons();
  
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
  Common::SafePtr<EulerSATerm> _model;

}; // end of class Euler2DSAPuvtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_Euler2DSAPuvtToCons_hh
