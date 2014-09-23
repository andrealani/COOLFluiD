#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobAnalytic_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobAnalytic_hh

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
class FVMCC_ComputeRhsJacobAnalytic : public FVMCC_ComputeRhsJacob {
public:
  
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobAnalytic(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCC_ComputeRhsJacobAnalytic();

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
  virtual void computeBothJacobTerms();
   
  /**
   * Compute only the jacobian contribution to one of the two states 
   * (this is used in parallel computing)
   */
  virtual void computeJacobTerm(const CFuint idx);

  /// axisymmetric jacobian term
  void axiJacobTerm(CFuint idx, RealMatrix& convJacobL, RealMatrix& convJacobR, CFreal factor);
  
  /// non axisymmetric jacobian term
  void noaxiJacobTerm(CFuint idx, RealMatrix& convJacobL, RealMatrix& convJacobR, CFreal factor);
  
  /// axisymmetric jacobian terms (left and right)
  void axiBothJacobTerms(RealMatrix& convJacobL, RealMatrix& convJacobR);
  
  /// non axisymmetric jacobian terms (left and right)
  void noaxiBothJacobTerms(RealMatrix& convJacobL, RealMatrix& convJacobR);
  
private:
  
  /// temporary left jacobian
  RealMatrix _tmpJacobL;
  
  /// temporary right jacobian
  RealMatrix _tmpJacobR;
    
}; // class FVMCC_ComputeRhsJacobAnalytic

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobAnalytic_hh
