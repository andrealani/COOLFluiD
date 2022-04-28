#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionMHD_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"
#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for MHD
   *
   * @author Ray Vandenhoeck
   *
   */
class DiffBndCorrectionsRHSFluxReconstructionMHD : public DiffBndCorrectionsRHSFluxReconstruction {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSFluxReconstructionMHD(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSFluxReconstructionMHD();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: // functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
//  /// compute the interface flux
//  void computeInterfaceFlxCorrection();
  
protected: //data
    
  /// matrix to store the state terms needed for the gradients (p, u, v, T)
  RealMatrix m_tempGradTerm;
  
  /// diffusive variable set
  Common::SafePtr< Physics::MHD::MHDProjectionDiffVarSet > m_diffusiveVarSet;
  
  /// element states of the neighbor in the correct format
  std::vector< RealVector* > m_tempStates;

}; // end of class DiffBndCorrectionsRHSFluxReconstructionMHD

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionMHD_hh
