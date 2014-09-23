#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRHS.hh"
#include "Common/CFMap.hh"
#include "Framework/BlockAccumulator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {
      class FVMCC_BC;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_ComputeRhsJacob : public FVMCC_ComputeRHS {
public:
    
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacob();

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
   * Compute only the jacobian contribution to one of the two states
   * (this is used in parallel computing)
   */
  virtual void computeJacobTerm(const CFuint idx);

protected:
  
  /**
   * Compute the jacobian contribution of the current (boundary) face
   */
  virtual void computeBoundaryJacobianTerm();
  
  /// Initialize the computation of RHS
  virtual void initializeComputationRHS();
  
  /// Compute the jacobian of the RHS
  virtual void computeRHSJacobian();
  
  /// Finalize the computation of RHS
  virtual void finalizeComputationRHS();
  
  /**
   * Compute the source term and its numerical jacobian
   * @param ist ID of the source term
   */
  void addSourceTermNumJacob(Framework::GeometricEntity* const cell, CFuint idx)
  {
    for (CFuint i = 0; i < _stNumJacobIDs.size(); ++i) {
      const CFuint ist = _stNumJacobIDs[i];
      (*_stComputers)[ist]->computeSource(cell, _pertSource[ist], _dummyJacob);
      // compute the finite difference derivative of the source term
      _numericalJacob->computeDerivative(_source[idx][ist], _pertSource[ist], _sourceDiff);
      _sourceDiffSum += _sourceDiff;
    }
  }
  
  /// Add jacobian term
  void addJacobTerm(CFuint idx, CFuint iVar, CFuint iCell, Framework::BlockAccumulator *const acc);
  
  /// Add both jacobian terms (left and right)
  void addBothJacobTerms(CFuint iVar, CFuint iCell);
  
  /// Add the analytical source term contribution
  void addAnalyticSourceTermJacob(CFuint idx, Framework::BlockAccumulator *const acc)
  {    
    for (CFuint i = 0; i < _stAnJacobIDs.size(); ++i) {
      const CFuint ist = _stAnJacobIDs[i];
      RealMatrix& sourceJacob = _sourceJacobian[idx][ist];
      sourceJacob *= _upStFactor[idx];
      acc->addValuesM(idx, idx, sourceJacob); 
      
      //  if (_currFace->getNeighborGeo(idx)->getState(0)->getLocalID() == 100) {
      // 	std::cout << sourceJacob << std::endl;
      //       }
    }
  }
    
  /// Compute convective and diffusive fluxes
  virtual void computeConvDiffFluxes(CFuint iVar, CFuint iCell);
  
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
    
  /// perturbed source terms
  std::vector<RealVector> _pertSource;
  
  /// source term delta
  RealVector _sourceDiff;
  
  /// sum of all source term deltas
  RealVector _sourceDiffSum;
  
  /// dummy jacobian matrix
  RealMatrix _dummyJacob;
  
}; // class FVMCC_ComputeRhsJacob

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobNumerics_hh
