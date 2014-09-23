#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQConsToSymm1InRef_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQConsToSymm1InRef_hh

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
 * This class represents a transformer of variables:
 * dQ/dU
 *
 * @author Andrea Lani
 * @modified Jesus Garicano Mena
 */
class Euler2DNEQConsToSymm1InRef : public Framework::VarSetMatrixTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQConsToSymm1InRef(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQConsToSymm1InRef();

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
  
}; // end of class Euler2DNEQConsToSymm1InRef

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQConsToSymm1InRef_hh
