#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for implicit schemes for NS
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   * @author Rayan Dhib
   *
   */
class DiffBndCorrectionsRHSJacobFluxReconstructionNS : public DiffBndCorrectionsRHSJacobFluxReconstruction {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSJacobFluxReconstructionNS(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSJacobFluxReconstructionNS();
  
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
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
    
}; // end of class DiffBndCorrectionsRHSJacobFluxReconstructionNS

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionNS_hh
