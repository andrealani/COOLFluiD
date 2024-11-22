#ifndef COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhovtXiToRhoivtTv_hh
#define COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhovtXiToRhoivtTv_hh

//////////////////////////////////////////////////////////////////////

#include "NEQ/Euler3DNEQRhovtXiToRhoivtTv.hh"

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
class PoissonNEQ3DRhovtXiToRhoivtTv : public NEQ::Euler3DNEQRhovtXiToRhoivtTv {
public:
  
  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> PTERM;
  
  /**
   * Default constructor without arguments
   */
  PoissonNEQ3DRhovtXiToRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~PoissonNEQ3DRhovtXiToRhoivtTv();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
    
}; // end of class PoissonNEQ3DRhovtXiToRhoivtTv

//////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_PoissonNEQ_PoissonNEQ3DRhovtXiToRhoivtTv_hh
