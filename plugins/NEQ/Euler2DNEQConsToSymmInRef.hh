#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQConsToSymmInRef_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQConsToSymmInRef_hh

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
class Euler2DNEQConsToSymmInRef : public Framework::VarSetMatrixTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQConsToSymmInRef(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQConsToSymmInRef();

  /**
   * Set the transformation matrix from reference values
   */
  void setMatrixFromRef();

private:
  //
  void computeDirect(RealMatrix& mat);
  void computeInverse(RealMatrix& mat);
  //
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
   
}; // end of class Euler2DNEQConsToSymmInRef

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQConsToSymmInRef_hh
