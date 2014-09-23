#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobCoupling_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobCoupling_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRHS.hh"
#include "Common/CFMap.hh"
#include "Framework/BlockAccumulator.hh"

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
class FVMCC_ComputeRhsJacobCoupling : public FVMCC_ComputeRHS {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobCoupling(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacobCoupling();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
protected:

  /// Initialize the computation of RHS
  virtual void initializeComputationRHS();
  
  /// Compute the jacobian of the RHS
  virtual void computeRHSJacobian();
  
  /// Finalize the computation of RHS
  virtual void finalizeComputationRHS();
  
  /**
   * Compute the jacobian contribution of the current (internal) face
   */
  virtual void computeJacobianTerm();
  
  /**
   * Compute the jacobian contribution of the current (internal) face
   * for the corresponding set of LSS equations
   */
  virtual void computeBothJacobTerms();
  
  /**
   * Compute only the jacobian contribution to one of the two states 
   * (this is used in parallel computing)
   */
  virtual void computeJacobTerm(CFuint idx);
  
  /**
   * Compute the jacobian contribution of the current (boundary) face
   */
  virtual void computeBoundaryJacobianTerm();
  
  /// Compute convective and diffusive fluxes
  virtual void computeConvDiffFluxes(CFuint iVar, CFuint iCell);
  
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
      _numericalJacob->computeDerivative(_source[idx][ist].slice(_start, _nbEqs), 
					 _pertSource[ist].slice(_start, _nbEqs),
					 _sourceDiff.slice(_start, _nbEqs));
      _sourceDiffSum.slice(_start, _nbEqs) += _sourceDiff.slice(_start, _nbEqs);
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
      RealMatrix& sourceJacob = _sourceJacobian[idx][_stAnJacobIDs[i]];
      _tmpAnalytMatrix[_iLSS].slice(0,0, _nbEqs, _nbEqs) = sourceJacob.slice(_start,_start, _nbEqs, _nbEqs);
      _tmpAnalytMatrix[_iLSS] *= _upStFactor[idx]; 
      acc->addValuesM(idx, idx, _tmpAnalytMatrix[_iLSS]);
    }
  }
  
  /// set the data corresponding to the current subsystem of equations
  void setCurrentSubsystemData(const std::vector<CFuint>& currEqs) 
  {
    using namespace COOLFluiD::Framework;
    
    // number of equations in the current subsystem
    _nbEqs = currEqs.size();
    // ID corresponding to the first entry in the current subsystem
    _start = currEqs[0];
    // set the equation subsystem descriptor
    PhysicalModelStack::getActive()->setEquationSubSysDescriptor(_start, _nbEqs, _iLSS);
    CFLog(DEBUG_MIN, "FVMCC_ComputeRhsJacobCoupling::setEquationSubSysDescriptor() => (" 
	  << _nbEqs << "," << _start << "," << _iLSS <<   ")\n");
  }
  
protected:
  
  /// pointer to the linear system solver
  Common::SafePtr<Framework::NumericalJacobian> _numericalJacob;
  
  /// current number of equations
  CFuint _nbEqs;
  
  /// ID of the current starting variables
  CFuint _start;
  
  /// ID of the current LinearSystemSolver
  CFuint _iLSS;
  
  /// vector for the temporary flux corresponding to a perturbed state
  RealVector _pertFlux;

  /// vector for the temporary flux finite difference
  RealVector _fluxDiff;
  
  /// vector for the original state values
  RealVector _origState;
  
  /// acquaintance of all the linear system solvers
  std::vector<Common::SafePtr<Framework::LinearSystemSolver> > _lss;
  
  // accumulator for all the LSSMatrix's
  std::vector<Framework::BlockAccumulator*> _acc;
      
  // accumulator for all the LSSMatrix's for boundary contribution
  std::vector<Framework::BlockAccumulator*> _bAcc;
  
  // array of the equation IDs for each LSS
  std::vector<Common::SafePtr<std::vector<CFuint> > > _equations;
  
  /// perturbed source terms
  std::vector<RealVector> _pertSource;
  
  /// source term delta
  RealVector _sourceDiff;
  
  /// sum of all source term deltas
  RealVector _sourceDiffSum;
  
  /// dummy jacobian matrix
  RealMatrix _dummyJacob;
  
  /// temporary analytical matrix
  std::vector<RealMatrix> _tmpAnalytMatrix;
  
}; // class FVMCC_ComputeRhsJacobCoupling
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume
    
  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobCoupling_hh
