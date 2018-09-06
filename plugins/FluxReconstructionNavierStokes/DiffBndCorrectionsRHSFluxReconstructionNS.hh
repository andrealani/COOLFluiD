#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for NS
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
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;
  
  /// element states of the neighbor in the correct format
  std::vector< RealVector* > m_tempStates;

}; // end of class DiffBndCorrectionsRHSFluxReconstructionNS

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionNS_hh
