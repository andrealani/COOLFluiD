#ifndef COOLFluiD_Physics_ICP_ICPNEQ2DPivtTvToRhoivtTv_hh
#define COOLFluiD_Physics_ICP_ICPNEQ2DPivtTvToRhoivtTv_hh

//////////////////////////////////////////////////////////////////////

#include "NEQ/Euler2DNEQPivtTvToRhoivtTv.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace ICP {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [p_i u v T Tv] to [rho_i u v T Tv] variables
 *
 * @author Andrea Lani
 */
class ICPNEQ2DPivtTvToRhoivtTv : public NEQ::Euler2DNEQPivtTvToRhoivtTv {
public:
  
  typedef Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> PTERM;
  
  /**
   * Default constructor without arguments
   */
  ICPNEQ2DPivtTvToRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~ICPNEQ2DPivtTvToRhoivtTv();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
    
}; // end of class ICPNEQ2DPivtTvToRhoivtTv

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPNEQ2DPivtTvToRhoivtTv_hh
