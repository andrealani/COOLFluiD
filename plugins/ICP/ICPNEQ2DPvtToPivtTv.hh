#ifndef COOLFluiD_Physics_ICP_ICPNEQ2DPvtToPivtTv_hh
#define COOLFluiD_Physics_ICP_ICPNEQ2DPvtToPivtTv_hh

//////////////////////////////////////////////////////////////////////

#include "NEQ/Euler2DNEQPuvtToPivtTv.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace ICP {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [p u v t] to [p_i u v T Tv] variables
 *
 * @author Andrea Lani
 */
class ICPNEQ2DPvtToPivtTv : public NEQ::Euler2DNEQPuvtToPivtTv {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  ICPNEQ2DPvtToPivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~ICPNEQ2DPvtToPivtTv();

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
    throw Common::NotImplementedException (FromHere(),"Euler2DNEQPvtToPivtTv::transformFromRef()");
  }
  
}; // end of class ICPNEQ2DPvtToPivtTv

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPNEQ2DPvtToPivtTv_hh
