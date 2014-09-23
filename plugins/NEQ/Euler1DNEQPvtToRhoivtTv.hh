#ifndef COOLFluiD_Physics_NEQ_Euler1DNEQPvtToRhoivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler1DNEQPvtToRhoivtTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from primitive
 * [p u T] to conservative variables with LTE and Demixing
 *
 * @author Alessandro Munafò
 *
 */
class Euler1DNEQPvtToRhoivtTv : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
    
  /**
   * Default constructor without arguments
   */
  Euler1DNEQPvtToRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler1DNEQPvtToRhoivtTv();

  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result)
  {
    throw Common::NotImplementedException (FromHere(),"Euler1DNEQPvtToRhoivtTv::transformFromRef()");
  }
    
private:

   /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
    
  /// Vector storing the elemental composition
  RealVector _ye;

  /// Vector storing the elemental composition
  RealVector _xe;

}; // end of class Euler1DNEQPvtToRhoivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler1DNEQPuvtToRhoivtTv_hh
