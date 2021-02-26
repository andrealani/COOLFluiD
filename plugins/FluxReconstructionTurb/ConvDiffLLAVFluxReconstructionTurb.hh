// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVFluxReconstructionTurb_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVFluxReconstructionTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/ConvDiffLLAVFluxReconstructionNS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes the contribution of the 
 * convective, diffusive and LLAV parts to the RHS for a Flux Reconstruction scheme for implicit 
 * time marching for RANS: same as NS, but added wall distance for diffusive
 * flux computations
 * 
 * @author Ray Vandenhoeck
 */
class ConvDiffLLAVFluxReconstructionTurb : public ConvDiffLLAVFluxReconstructionNS {

public: // functions

  /// Constructor
  explicit ConvDiffLLAVFluxReconstructionTurb(const std::string& name);

  /// Destructor
  virtual ~ConvDiffLLAVFluxReconstructionTurb() {}
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected: //functions
    
    /**
   * compute the unperturbed cell diffusive residuals
   * @pre m_faceTermComputers->computeDiffFaceTermAndUpdateCoefContributions
   * @pre setCellsData()
   */
  virtual void computeUnpertCellDiffResiduals(const CFuint side);
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
protected: // data

  /// handle to the wall distance
  Framework::DataSocketSink<CFreal> socket_wallDistance;
  
  /// idx of closest sol to each flx
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlxIdx;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  /// Physical data temporary vector
  RealVector m_pData2;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVFluxReconstructionTurb_hh
