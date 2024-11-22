#ifndef COOLFluiD_Physics_NEQ_Euler3DNEQRhovtXiToRhoivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler3DNEQRhovtXiToRhoivtTv_hh

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
 * [rho u v w T] to conservative variables with LTE and Demixing
 *
 * @author Vatsalya Sharma
 *
 */
class Euler3DNEQRhovtXiToRhoivtTv : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
    
  /**
   * Default constructor without arguments
   */
  Euler3DNEQRhovtXiToRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler3DNEQRhovtXiToRhoivtTv();

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
    throw Common::NotImplementedException (FromHere(),"Euler3DNEQRhovtXiToRhoivtTv::transformFromRef()");
  }
    
protected: //Vatsalya: was private earlier.

  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
  
  /// Vector storing the mass fractions
  RealVector _ye;
  
  /// Vector storing the molar fractions
  RealVector _xe;

}; // end of class Euler1DNEQPvtToRhoivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler3DNEQRhovtXiToRhoivtTv_hh
