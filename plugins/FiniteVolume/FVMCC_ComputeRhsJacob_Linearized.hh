#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_Linearized_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_Linearized_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace Numerics {

    namespace FiniteVolume {
      class PolyReconstructorLin;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_ComputeRhsJacob_Linearized : public FVMCC_ComputeRhsJacob {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacob_Linearized(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacob_Linearized();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Execute Processing actions
   */
  virtual void execute();

protected:

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

  /**
   * Compute the jacobian contribution of the current (boundary) face
   */
  virtual void computeBoundaryJacobianTerm();

  /**
   * Set the face integrator data
   */
  virtual void setFaceIntegratorData();

  /**
   * Set the face normal
   */
  virtual void setFaceNormal();

  /**
   * Compute the physical data in the states
   */
  void computeStatesData2()
  {
    _reconstrVar->computeStatesData
      (_rStates, _rExtraVars, 4*_nbQPointsInFace[_faceIdx]);
  }

  /**
   * Back up the physical data in the states
   */
  void computeAndBackUpStatesData2()
  {
    _reconstrVar->computeAndBackUpStatesData
      (_rStates, _rExtraVars, 4*_nbQPointsInFace[_faceIdx]);
  }

  /**
   * Restore the backed up left states
   */
  void restoreLeftStates2()
  {
    _reconstrVar->restoreBackupStates
      (_fluxData2->leftValues, 2*_nbQPointsInFace[_faceIdx]);
  }

  /**
   * Restore the backed up right states
   */
  void restoreRightStates2()
  {
    _reconstrVar->restoreBackupStates
      (_fluxData2->rightValues, 2*_nbQPointsInFace[_faceIdx]);
  }

  /**
   * Compute the left perturbed states data
   */
  void computePerturbedLeftStatesData2(const CFuint iVar)
  {
    _reconstrVar->computePerturbedStatesData
      (_fluxData2->leftValues, 2*_nbQPointsInFace[_faceIdx], iVar);
  }

  /**
   * Compute the right perturbed states data
   */
  void computePerturbedRightStatesData2(const CFuint iVar)
  {
    _reconstrVar->computePerturbedStatesData
      (_fluxData2->rightValues, 2*_nbQPointsInFace[_faceIdx], iVar);
  }

protected:

  /// storage for temporary extrapolated solution
  /// in the quadrature points
  Common::SafePtr<Framework::FluxSplitterData> _fluxData2;
  
  /// specific polynomial reconstructor
  Common::SafePtr<PolyReconstructorLin> _polyRecLin;  
  
}; // class FVMCC_ComputeRhsJacob_Linearized

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_LinearizedNumerics_hh
