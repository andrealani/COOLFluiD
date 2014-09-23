#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_CrankNichLimComputeRhs_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_CrankNichLimComputeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_ComputeRhsJacob.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Thomas Wuilbaut
 *
 */
class FVMCC_CrankNichLimComputeRhs : public FVMCC_ComputeRhsJacob {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_CrankNichLimComputeRhs(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_CrankNichLimComputeRhs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

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

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Un Setup private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * Get the factor multiplying the residual
   */
  CFreal getResFactorCN(const CFuint iEq) const
  {
    Framework::DataHandle<CFreal> timeLimiter = socket_timeLimiter.getDataHandle();
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    const CFuint state0ID = _currFace->getState(0)->getLocalID();
    const CFuint state1ID = _currFace->getState(1)->getLocalID();
    const CFreal timeLimValue  = 0.5*(timeLimiter(state0ID, iEq, nbEqs) + timeLimiter(state1ID, iEq, nbEqs));
    return (1. - timeLimValue*0.5);
  }
  
  /// Compute the update factors for left and right states in non axisymmetric case
  void computeNoAxiUpFactorsVec()
  {   
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    for(CFuint iEq=0; iEq<nbEqs; ++iEq){
      _upFactorVec[0][iEq] = getResFactorCN(iEq);
      _upStFactorVec[0][iEq] = _upStFactorVec[1][iEq] = -getResFactorCN(iEq);
    }
    
    _upFactorVec[1] = -1.0;
  }
  
  /// Compute the update factors for left and right states in axisymmetric case
  void computeAxiUpFactorsVec()
  { 
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    for(CFuint iEq=0; iEq<nbEqs; ++iEq){
      const CFreal lim = getResFactorCN(iEq);
      _upFactorVec[0][iEq] = _rMid*_invr[0]*lim;
      _upStFactorVec[0][iEq] = -lim*_invr[0];
      _upStFactorVec[1][iEq] = -lim*_invr[1];
    }
    
    _upFactorVec[1] = (-_invr[1]/_invr[0]);
  }
  
  /// Add the analytical source term contribution
  void addAnalyticSourceTermJacobCN(CFuint idx, Framework::BlockAccumulator *const acc)
  {    
   //  for (CFuint i = 0; i < _stAnJacobIDs.size(); ++i) {
//       const CFuint ist = _stAnJacobIDs[i];
//       RealMatrix& sourceJacob = _sourceJacobian[idx][ist];
//       sourceJacob *= _upStFactorVec[idx]; //this is a matrix*vector !!
//       acc->addValues(idx, idx, sourceJacob); 
//     }
    
    std::cout << "FVMCC_CrankNichLimComputeRhs::addAnalyticSourceTermJacobCN() to be fixed" << std::endl; abort();
  }
   
  /// Add both jacobian terms (left and right)
  void addBothJacobTermsCN(CFuint iVar, CFuint iCell);
  
  /// Add jacobian term
  void addJacobTermCN(CFuint idx, CFuint iVar, CFuint iCell, 
		      Framework::BlockAccumulator *const acc);
  
  /**
   * Compute the contribution of the current face to the RHS
   */
  virtual void updateRHS();
  
protected:
  
  /// storage of the (time) limiter
  Framework::DataSocketSink< CFreal> socket_timeLimiter;
  
  /// update factor for fluxes LR
  std::vector<RealVector> _upFactorVec;
  
  /// update factor for source terms LR 
  std::vector<RealVector> _upStFactorVec;
  
}; // class FVMCC_CrankNichLimComputeRhs
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_CrankNichLimComputeRhsNumerics_hh
