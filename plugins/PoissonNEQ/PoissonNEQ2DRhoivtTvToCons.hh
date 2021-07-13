#ifndef COOLFluiD_Physics_PoissonNEQ_PoissonNEQ2DRhoivtTvToCons_hh
#define COOLFluiD_Physics_PoissonNEQ_PoissonNEQ2DRhoivtTvToCons_hh

//////////////////////////////////////////////////////////////////////

#include "NEQ/Euler2DNEQRhoivtTvToCons.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [rho_i u v T Tv] to conservative variables
 *
 * @author Andrea Lani
 */
class PoissonNEQ2DRhoivtTvToCons : public NEQ::Euler2DNEQRhoivtTvToCons {
public:
  
  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> PTERM;
  
  /**
   * Default constructor without arguments
   */
  PoissonNEQ2DRhoivtTvToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~PoissonNEQ2DRhoivtTvToCons();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
    
}; // end of class PoissonNEQ2DRhoivtTvToCons

//////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_PoissonNEQ_PoissonNEQ2DRhoivtTvToCons_hh
