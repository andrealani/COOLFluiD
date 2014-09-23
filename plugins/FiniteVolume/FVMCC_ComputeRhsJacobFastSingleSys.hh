#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFastSingleSys_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFastSingleSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRhsJacobFastCoupling.hh"

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
class FVMCC_ComputeRhsJacobFastSingleSys : public FVMCC_ComputeRhsJacobFastCoupling {
public:
  
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobFastSingleSys(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacobFastSingleSys();
  
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
  
}; // class FVMCC_ComputeRhsJacobFastSingleSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFastSingleSys_hh
