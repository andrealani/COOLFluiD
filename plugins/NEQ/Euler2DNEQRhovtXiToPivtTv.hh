#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQRhovtXiToPivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQRhovtXiToPivtTv_hh

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
 * This class represents a transformer of variables from [rho v t X_i] to 
 * [p_i v t Tv] variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DNEQRhovtXiToPivtTv : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
    
  /**
   * Default constructor without arguments
   */
  Euler2DNEQRhovtXiToPivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQRhovtXiToPivtTv();

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
    throw Common::NotImplementedException (FromHere(),"Euler2DNEQRhovtXiToPivtTv::transformFromRef()");
  }
    
private:

  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
  
  /// Vector storing the mass fractions
  RealVector _ye;
  
  /// Vector storing the molar fractions
  RealVector _xe;
  
  /// Vector storing the temperatures
  RealVector _tvDim;
  
  /// array with all species molar masses
  RealVector _Rspecies;
  
}; // end of class Euler1DNEQPvtToPivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQRhovtXiToPivtTv_hh
