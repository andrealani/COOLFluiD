#ifndef COOLFluiD_Physics_ICP_ICPNEQ2DPuvtToRhoivtTv_hh
#define COOLFluiD_Physics_ICP_ICPNEQ2DPuvtToRhoivtTv_hh

//////////////////////////////////////////////////////////////////////

#include "NEQ/Euler2DNEQPuvtToRhoivtTv.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace ICP {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [p u v t] [rho_i u v T Tv] to conservative variables
 *
 * @author Andrea Lani
 */
class ICPNEQ2DPuvtToRhoivtTv : public NEQ::Euler2DNEQPuvtToRhoivtTv {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  ICPNEQ2DPuvtToRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~ICPNEQ2DPuvtToRhoivtTv();

  /**
   * Transform a set of state vectors into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
 /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DNEQPuvtToRhoivtTv::transformFromRef()");
  }
  
}; // end of class ICPNEQ2DPuvtToRhoivtTv

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPNEQ2DPuvtToRhoivtTv_hh
