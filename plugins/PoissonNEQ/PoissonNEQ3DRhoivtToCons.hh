#ifndef COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhoivtToCons_hh
#define COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhoivtToCons_hh

//////////////////////////////////////////////////////////////////////

#include "NEQ/Euler3DNEQRhoivtToCons.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [rho_i u v T ] to conservative variables
 *
 * @author Andrea Lani
 */
class PoissonNEQ3DRhoivtToCons : public NEQ::Euler3DNEQRhoivtToCons {
public:
  
  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> PTERM;
  
  /**
   * Default constructor without arguments
   */
  PoissonNEQ3DRhoivtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~PoissonNEQ3DRhoivtToCons();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
    
}; // end of class PoissonNEQ3DRhoivtToCons

//////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhoivtToCons_hh
