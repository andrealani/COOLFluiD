#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionTurb_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSFluxReconstructionNS.hh"

#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for diffusive terms to the RHS for implicit schemes for Turb
   *
   * @author Ray Vandenhoeck
   *
   */
class DiffBndCorrectionsRHSFluxReconstructionTurb : public DiffBndCorrectionsRHSFluxReconstructionNS {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSFluxReconstructionTurb(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSFluxReconstructionTurb();
  
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
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected: // functions

  
protected: // data
    
    /// handle to the wall distance
  Framework::DataSocketSink<CFreal> socket_wallDistance;
  
  /// idx of closest sol to each flx
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlxIdx;
  
  Common::SafePtr< Physics::KOmega::NavierStokes2DKLogOmega > m_navierStokesVarSetTurb;
  
  Common::SafePtr< Physics::KOmega::NavierStokes3DKLogOmega > m_navierStokesVarSetTurb3D;

    
}; // end of class DiffBndCorrectionsRHSFluxReconstructionTurb

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionTurb_hh
