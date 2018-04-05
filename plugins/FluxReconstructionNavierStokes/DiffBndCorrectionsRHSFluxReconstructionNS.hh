#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for convective terms to the RHS
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   *
   */
class DiffBndCorrectionsRHSFluxReconstructionNS : public DiffBndCorrectionsRHSFluxReconstruction {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSFluxReconstructionNS(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSFluxReconstructionNS();

protected: // functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();

}; // end of class DiffBndCorrectionsRHSFluxReconstructionNS

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNS_hh
