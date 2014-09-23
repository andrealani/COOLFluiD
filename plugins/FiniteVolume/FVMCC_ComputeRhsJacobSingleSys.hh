#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobSingleSys_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobSingleSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRhsJacobCoupling.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {
  
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_ComputeRhsJacobSingleSys : public FVMCC_ComputeRhsJacobCoupling {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobSingleSys(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacobSingleSys();
  
  /**
   * Execute Processing actions
   */
  virtual void execute();
  
protected:
  
  /**
   * Compute the jacobian contribution of the current (internal) face
   */
  virtual void computeJacobianTerm();
  
  /**
   * Compute the source term contribution
   */
  virtual void computeSourceTermContribution();
    
protected:
  
  /// current LSS idx
  CFuint _currLSSIdx;
  
}; // class FVMCC_ComputeRhsJacobSingleSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobSingleSys_hh
