#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToConsInRef_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToConsInRef_hh

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
class Euler2DNEQRhoivtTvToConsInRef : public Framework::VarSetMatrixTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQRhoivtTvToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQRhoivtTvToConsInRef();

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
  
}; // end of class Euler2DNEQRhoivtTvToConsInRef

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQRhoivtTvToConsInRef_hh
