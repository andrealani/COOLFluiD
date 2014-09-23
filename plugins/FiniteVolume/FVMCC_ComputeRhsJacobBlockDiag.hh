#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobBlockDiag_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobBlockDiag_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRHS.hh"
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
 * @author Andrea Lani
 * @author Jiri Simonek
 *
 */
class FVMCC_ComputeRhsJacobBlockDiag : public FVMCC_ComputeRHS {
public:
    
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobBlockDiag(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacobBlockDiag();

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
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

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
   * Compute the jacobian contribution of the current (internal) face
   */
  virtual void computeBothJacobTerms();

  /**
   * Compute the source term contribution
   */
  virtual void computeSourceTermContribution();

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

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> _lss;

  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;
  
  /// vector for the temporary flux corresponding to a perturbed state
  RealVector _pertFlux;

  /// vector for the temporary flux finite difference
  RealVector _fluxDiff;

  /// vector for the original state values
  RealVector _origState;

  // accumulator for LSSMatrix
  std::auto_ptr<Framework::BlockAccumulator> _acc;

  // accumulator for LSSMatrix for boundary contribution
  std::auto_ptr<Framework::BlockAccumulator> _bAcc;

  /// array storing temporary Element nodes.
  std::vector<Framework::Node*> _cellNodes;

  /// axisymmetric flux
  RealVector _axiFlux;
    
}; // class FVMCC_ComputeRhsJacobBlockDiag

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobBlockDiagNumerics_hh
