#ifndef COOLFluiD_Physics_SA_EulerSAPvtToCons_hh
#define COOLFluiD_Physics_SA_EulerSAPvtToCons_hh

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
 * [p u v w T Nuitil] to conservative variables [rho rhoU rhoV rhoW rhoE rhoNuitil]
 *
 * @author Joao Pinto
 * @author Thomas Wuilbaut
 * @modified Christos Gkoudesnes
 * @modified Andrea Lani
 *
 */
template <typename BASE>
class EulerSAPvtToCons : public BASE {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> EulerSATerm;
  
  /**
   * Constructor
   */
  EulerSAPvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);
  
  /**
   * Default destructor
   */
  virtual ~EulerSAPvtToCons();
  
  /**
   * Transform a state into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
  
private:
  
  /// acquaintance of the PhysicalModel
  Common::SafePtr<EulerSATerm> m_modelSA;
  
}; // end of class EulerSAPvtToCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "SA/EulerSAPvtToCons.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_SA_EulerSAPvtToCons_hh
