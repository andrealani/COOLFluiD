#ifndef COOLFluiD_Physics_ICP_ICPNEQ2DPvtToPivtTv_hh
#define COOLFluiD_Physics_ICP_ICPNEQ2DPvtToPivtTv_hh

//////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace ICP {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from
 * [p_i u v T Tv] to conservative variables
 *
 * @author Andrea Lani
 */
class ICPNEQ2DPvtToPivtTv : public Framework::VarSetTransformer {
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
  void transformFromRef(const RealVector& data, Framework::State& result)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DNEQPuvtToRhoivtTv::transformFromRef()");
  }
           
private:
  
  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
    
  /// Vector storing the elemental composition
  RealVector _ye;

  /// Vector storing the elemental composition
  RealVector _xe;
  
  /// Vector storing the species masses
  RealVector _masses;
  
}; // end of class ICPNEQ2DPvtToPivtTv

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPNEQ2DPvtToPivtTv_hh
