#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQPivtTvToRhoivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQPivtTvToRhoivtTv_hh

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
 * This class represents a transformer of variables from [pi u v T Tvi] to 
 * [rhoi u v T Tvi]
 *
 * @author Andrea Lani
 *
 */
class Euler2DNEQPivtTvToRhoivtTv : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQPivtTvToRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~Euler2DNEQPivtTvToRhoivtTv();
  
  /**
   * Transform a state into another one
   */
  virtual void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  virtual void transformFromRef(const RealVector& data, Framework::State& result);
  
protected:
  
  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
  
}; // end of class Euler2DNEQPivtTvToRhoivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQPivtTvToRhoivtTv_hh
