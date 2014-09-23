#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFastCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFastCoupling_hh

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
class FVMCC_ComputeRhsJacobFastCoupling : public FVMCC_ComputeRhsJacobCoupling {
public:
  
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobFastCoupling(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCC_ComputeRhsJacobFastCoupling();

  /**
   * Execute Processing actions
   */
  virtual void execute();
  
protected:
  
  /**
   * Compute the jacobian contribution of the current (internal) face
   * for the corresponding set of LSS equations
   */
  virtual void computeBothJacobTerms(CFuint iLSS);
  
  /**
   * Compute only the jacobian contribution to one of the two states 
   * (this is used in parallel computing)
   */
  virtual void computeJacobTerm(CFuint iLSS, CFuint idx);
  
  /**
   * Compute the jacobian contribution of the current (boundary) face for 
   * the current subsystem of equations
   */
  virtual void computeBoundaryJacobianTerm(CFuint iLSS);
  
}; // class FVMCC_ComputeRhsJacobFastCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFastCoupling_hh
