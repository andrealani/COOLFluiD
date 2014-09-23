#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsLinearizedJacob_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsLinearizedJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRhsJacob_Linearized.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace Numerics {

    namespace FiniteVolume {
      class FVMCC_BC;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Thomas Wuilbaut
 *
 */
class FVMCC_ComputeRhsLinearizedJacob : public FVMCC_ComputeRhsJacob_Linearized {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsLinearizedJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsLinearizedJacob();


protected:

  /**
   * Compute the jacobian contribution of the current (internal) face
   */
  virtual void computeBothJacobTerms();

  /**
   * Compute only the jacobian contribution to one of the two states
   * (this is used in parallel computing)
   */
  virtual void computeJacobTerm(const CFuint idx);

protected:

  /**
   * Compute the jacobian contribution of the current (boundary) face
   */
  virtual void computeBoundaryJacobianTerm();


}; // class FVMCC_ComputeRhsLinearizedJacob

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsLinearizedJacobNumerics_hh
