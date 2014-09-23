#ifndef COOLFluiD_Numerics_FiniteVolume_NullComputeSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_NullComputeSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
    class ConvectiveVarSet;
    class DiffusiveVarSet;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

 /**
  * This class offers and abstract interface to computes the source
  * term for cell center FVM schemes.
  *
  * @author Andrea Lani
  *
  */
class NullComputeSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Constructor
   */
  NullComputeSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~NullComputeSourceTerm();
    
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
		     RealVector& source,
		     RealMatrix& jacobian);

  
}; // end of class NullComputeSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NullComputeSourceTerm_hh
