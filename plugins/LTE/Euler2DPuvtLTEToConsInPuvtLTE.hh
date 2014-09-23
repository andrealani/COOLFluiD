#ifndef COOLFluiD_Physics_LTE_Euler2DPuvtLTEToConsInPuvtLTE_hh
#define COOLFluiD_Physics_LTE_Euler2DPuvtLTEToConsInPuvtLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetMatrixTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class EulerTerm;
    }
    
    namespace LTE {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from conservative
 * to primitive [p u v T] starting from primitive variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DPuvtLTEToConsInPuvtLTE : public
Framework::VarSetMatrixTransformer {

public:

  /**
   * Default constructor without arguments
   */
  Euler2DPuvtLTEToConsInPuvtLTE
  (Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  virtual ~Euler2DPuvtLTEToConsInPuvtLTE();

  /**
   * Set the transformation matrix from a given state
   */
  virtual void setMatrix(const RealVector& state);

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

  /// acquaintance of the PhysicalModel
  Common::SafePtr<NavierStokes::EulerTerm> _model;

  /// tolerance for numerical derivation
  CFreal _tol;
  
  /// array to store density, enthalpy and energy
  RealVector _dhe;
  
}; // end of class Euler2DPuvtLTEToConsInPuvtLTE

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_Euler2DPuvtLTEToConsInPuvtLTE_hh
