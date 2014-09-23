#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQSymm1ToConsInRef_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQSymm1ToConsInRef_hh

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
 * This class implements the dU/dQ matrix
 * 
 *
 * @author Andrea Lani
 * @modified  Jesus Garicano Mena
 */
class Euler2DNEQSymm1ToConsInRef : public Framework::VarSetMatrixTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQSymm1ToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQSymm1ToConsInRef();

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
 
  /// array to store the mass fractions of elements
  RealVector _ys;
 
  /// useful coefficients
  RealVector _alpha;
  
  /// useful coefficients
  RealVector _RiGas;
  
  /// species molar masses
  RealVector _mmasses;
   
  /// f_i coefficients
  RealVector _fcoeff;
}; // end of class Euler2DNEQSymm1ToConsInRef

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQSymm1ToConsInRef_hh
