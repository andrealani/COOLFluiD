#ifndef COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhoivtTvToCons_hh
#define COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhoivtTvToCons_hh

//////////////////////////////////////////////////////////////////////

#include "NEQ/Euler3DNEQRhoivtTvToCons.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [rho_i u v w T Tv] to conservative variables
 *
 * @author VS
 */
class PoissonNEQ3DRhoivtTvToCons : public NEQ::Euler3DNEQRhoivtTvToCons {
public:
  
  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> PTERM;
  
  /**
   * Default constructor without arguments
   */
  PoissonNEQ3DRhoivtTvToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~PoissonNEQ3DRhoivtTvToCons();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
    
}; // end of class PoissonNEQ3DRhoivtTvToCons

//////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhoivtTvToCons_hh
