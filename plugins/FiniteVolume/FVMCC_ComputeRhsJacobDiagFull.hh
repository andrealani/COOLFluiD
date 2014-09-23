#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobDiagFull_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobDiagFull_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_ComputeRhsJacobDiag.hh"

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
class FVMCC_ComputeRhsJacobDiagFull : public FVMCC_ComputeRhsJacobDiag {
public:
  
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobDiagFull(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCC_ComputeRhsJacobDiagFull();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

private:
  
  /**
   * Compute the jacobian contribution of the current (internal) face
   */
  virtual void computeBothJacobTerms();
  
  /**
   * Compute only the jacobian contribution to one of the two states 
   * (this is used in parallel computing)
   */
  virtual void computeJacobTerm(const CFuint idx);
  
  /**
   * Compute the jacobian contribution of the current (boundary) face
   */
  virtual void computeBoundaryJacobianTerm();
  
private:
  
  /// temporary left jacobian
  RealMatrix _tmpJacob;
  
  /// average state
  Framework::State* _avState;
  
}; // class FVMCC_ComputeRhsJacobDiagFull

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobDiagFull_hh
