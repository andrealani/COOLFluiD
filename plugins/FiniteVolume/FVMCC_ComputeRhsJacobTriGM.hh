#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobTriGM_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobTriGM_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"
//#include "Common/CFMap.hh"

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
 * @author Jiri Simonek
 *
 */
class FVMCC_ComputeRhsJacobTriGM : public FVMCC_ComputeRhsJacob {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobTriGM(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacobTriGM();

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

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();

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
  
  /**
		 * Compute the source term and its analytical jacobian
		 * @param ist ID of the source term
	 */
		virtual void computeSourceTermAnalytJacob(CFuint ist);
  
  /**
		 * Compute the source term and its numerical jacobian
		 * @param ist ID of the source term
	 */
		virtual void computeSourceTermNumJacob(CFuint ist);	

	protected:

  /// storage of diagonal block matrices (jacobians) in block tridiagonal matrix
		Framework::DataSocketSource<CFreal> socket_diagMatrices;

  /// storage of under diagonal block matrices (jacobians) in block tridiagonal matrix
		Framework::DataSocketSource<CFreal> socket_underDiagMatrices;

  /// sotrage of above diagonal block matrices (jacobians) in block tridiagonal matrix
		Framework::DataSocketSource<CFreal> socket_aboveDiagMatrices;

  /// storage of the local updatable IDs or -1 (ghost) for all local states
		Framework::DataSocketSource<CFint> socket_upLocalIDsAll;

  /// matrix iterator for diagonal matrices
		RealMatrix _matIterDiag0;

  /// matrix iterator for diagonal matrices
		RealMatrix _matIterDiag1;

  /// matrix iterator for under-diagonal matrices
		RealMatrix _matIterUnderDiag0;

  /// matrix iterator for under-diagonal matrices
		RealMatrix _matIterUnderDiag1;

  /// matrix iterator for above-diagonal matrices
		RealMatrix _matIterAboveDiag0;

  /// matrix iterator for above-diagonal matrices
		RealMatrix _matIterAboveDiag1;

}; // class FVMCC_ComputeRhsJacobTriGM

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobTriGMNumerics_hh
