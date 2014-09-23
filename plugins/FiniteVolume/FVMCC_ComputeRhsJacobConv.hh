#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobConv_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobConv_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"

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
class FVMCC_ComputeRhsJacobConv : public FVMCC_ComputeRhsJacob {
public:
  
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobConv(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCC_ComputeRhsJacobConv();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

private:
  
  /**
   * Compute the jacobian contribution of the current (internal) face
   */
  void computeBothJacobTerms();
   
  /**
   * Compute only the jacobian contribution to one of the two states 
   * (this is used in parallel computing)
   */
  void computeJacobTerm(const CFuint idx);
  
private:
  
  /// temporary left jacobian
  RealMatrix _tmpJacobL;
  
  /// temporary right jacobian
  RealMatrix _tmpJacobR;
    
}; // class FVMCC_ComputeRhsJacobConv

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobConv_hh
