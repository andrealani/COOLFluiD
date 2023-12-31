#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionPoisson_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
    namespace Physics {
    namespace Poisson {
      class PoissonDiffVarSet;
      class PoissonConvVarSet;
    }
  }

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for implicit schemes for Poisson
   *
   * @author Ray Vandenhoeck
   *
   */
class DiffBndCorrectionsRHSJacobFluxReconstructionPoisson : public DiffBndCorrectionsRHSJacobFluxReconstruction {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSJacobFluxReconstructionPoisson(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSJacobFluxReconstructionPoisson();
  
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
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

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
  
    Common::SafePtr<Physics::Poisson::PoissonDiffVarSet> m_diffVarSetPoisson;
    Common::SafePtr<Physics::Poisson::PoissonConvVarSet> m_convVarSetPoisson;
    
}; // end of class DiffBndCorrectionsRHSJacobFluxReconstructionPoisson

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSJacobFluxReconstructionPoisson_hh
