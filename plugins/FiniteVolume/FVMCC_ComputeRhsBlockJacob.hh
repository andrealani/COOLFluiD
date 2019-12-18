#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsBlockJacob_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsBlockJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRHS.hh"
#include "Common/CFMap.hh"
#include "Framework/BlockAccumulatorBase.hh"

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
class FVMCC_ComputeRhsBlockJacob : public FVMCC_ComputeRHS {
public:
    
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsBlockJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsBlockJacob();

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
  void addJacobTerm(CFuint idx, CFuint idxUpFactor,
		    CFuint iVar, CFuint iCell, Framework::BlockAccumulatorBase *const acc);
  
  /// Add the analytical source term contribution
  void addAnalyticSourceTermJacob(CFuint idx, Framework::BlockAccumulatorBase *const acc)
  {    
    for (CFuint i = 0; i < _stAnJacobIDs.size(); ++i) {
      const CFuint ist = _stAnJacobIDs[i];
      RealMatrix& sourceJacob = _sourceJacobian[idx][ist];
      sourceJacob *= _upStFactor[idx];
      acc->addValuesM(idx, idx, sourceJacob); 
    }
  }
    
  /// Compute convective and diffusive fluxes
  virtual void computeConvDiffFluxes(CFuint iVar, CFuint iCell);

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
protected:
  
  /// storage of diagonal block matrix
  Framework::DataSocketSource<CFreal> socket_diagMatrix;
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;
  
  /// block accumulator
  std::auto_ptr<Framework::BlockAccumulatorBase> _acc;
  
  /// vector for the temporary flux corresponding to a perturbed state
  RealVector _pertFlux;
  
  /// vector for the temporary flux finite difference
  RealVector _fluxDiff;
  
  /// vector for the original state values
  RealVector _origState;
      
  /// perturbed source terms
  std::vector<RealVector> _pertSource;
  
  /// source term delta
  RealVector _sourceDiff;
  
  /// sum of all source term deltas
  RealVector _sourceDiffSum;
  
  /// dummy jacobian matrix
  RealMatrix _dummyJacob;
  
}; // class FVMCC_ComputeRhsBlockJacob

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsBlockJacob_hh
