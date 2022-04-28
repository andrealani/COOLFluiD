#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionMHD_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSJacobFluxReconstruction.hh"

#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for implicit schemes for MHD
   *
   * @author Ray Vandenhoeck   *
   */
class DiffBndCorrectionsRHSJacobFluxReconstructionMHD : public DiffBndCorrectionsRHSJacobFluxReconstruction {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSJacobFluxReconstructionMHD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSJacobFluxReconstructionMHD();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * unset up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void unsetup();

protected: // functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);
  
  /**
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);
  
  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
protected: // data

    // vector for temporary storing the states of flx pnts
    std::vector< std::vector< RealVector* > > m_tempStates;
    
    // vector for temporary storing the states of sol pnts
    std::vector< RealVector* > m_tempStatesSol;
    
    /// diffusive variable set
  Common::SafePtr< Physics::MHD::MHDProjectionDiffVarSet > m_diffusiveVarSet;
    
}; // end of class DiffBndCorrectionsRHSJacobFluxReconstructionMHD

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionMHD_hh
