#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQRoeVinokurToConsInRef_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQRoeVinokurToConsInRef_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to Puvt to consistent conservative variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DNEQRoeVinokurToConsInRef : public Framework::VarSetMatrixTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQRoeVinokurToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQRoeVinokurToConsInRef();

  /**
   * Set the transformation matrix from reference values
   */
  void setMatrixFromRef();

private:

  /**
   * Set the flag telling if the transformation is an identity one
   * @pre this method must be called during set up
   */
  bool getIsIdentityTransformation() const
  {
    return false;
  }

private: //data
  
  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;

  /// constant of each gas
  RealVector _RiGas;

  /// Holds dpdRhoi at bar conditions:
  RealVector _alpha_bar;

  /// Holds kappa = dpdRhoE at bar conditions:
  CFreal _beta_bar;

}; // end of class Euler2DNEQRoeVinokurToConsInRef

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQRoeVinokurToConsInRef_hh
