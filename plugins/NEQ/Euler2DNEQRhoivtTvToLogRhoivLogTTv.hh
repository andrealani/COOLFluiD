#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToLogRhoivLogTTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToLogRhoivLogTTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * Transformer from (rho_i, u, v, T, Tv) to (ln rho_i, u, v, ln T, ln Tv),
 * used by boundary conditions given in RhoivtTv variables.
 *
 * @author Rayan Dhib
 */
class Euler2DNEQRhoivtTvToLogRhoivLogTTv : public Framework::VarSetTransformer {
public:

  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;

  /// Constructor
  Euler2DNEQRhoivtTvToLogRhoivLogTTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /// Default destructor
  virtual ~Euler2DNEQRhoivtTvToLogRhoivLogTTv();

  /// Transform a state into another one
  virtual void transform(const Framework::State& state, Framework::State& result);

  /// Transform from physical data: not needed
  virtual void transformFromRef(const RealVector& data, Framework::State& result)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DNEQRhoivtTvToLogRhoivLogTTv::transformFromRef()");
  }

private:

  /// acquaintance of the model
  Common::SafePtr<NEQTerm> m_model;

}; // end of class Euler2DNEQRhoivtTvToLogRhoivLogTTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToLogRhoivLogTTv_hh
