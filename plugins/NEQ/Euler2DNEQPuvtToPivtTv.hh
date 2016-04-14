#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQPuvtToPivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQPuvtToPivtTv_hh

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
 * [p u v T] to conservative variables with LTE and Demixing
 *
 * @author Andrea Lani
 *
 */
class Euler2DNEQPuvtToPivtTv : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
    
  /**
   * Default constructor without arguments
   */
  Euler2DNEQPuvtToPivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQPuvtToPivtTv();

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
    throw Common::NotImplementedException (FromHere(),"Euler2DNEQPuvtToPivtTv::transformFromRef()");
  }
    
private:

   /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
    
  /// Vector storing the elemental composition
  RealVector _ye;

  /// Vector storing the elemental composition
  RealVector _xe;
  
  /// molar masses
  RealVector _masses;
  
}; // end of class Euler1DNEQPvtToPivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQPuvtToPivtTv_hh
