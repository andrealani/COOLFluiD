// Copyright (C) 2019 KU Leuven, von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionTurb_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/ConvDiffLLAVJacobFluxReconstructionNS.hh"

#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"


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
class ConvDiffLLAVJacobFluxReconstructionTurb : public ConvDiffLLAVJacobFluxReconstructionNS {

public: // functions

  /// Constructor
  explicit ConvDiffLLAVJacobFluxReconstructionTurb(const std::string& name);

  /// Destructor
  virtual ~ConvDiffLLAVJacobFluxReconstructionTurb() {}
  
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
  
  /// initialize the data needed for the jacobian
  virtual void initJacobianComputation();
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /// compute the Riemann flux jacobian numerically
  virtual void computeRiemannFluxJacobianNum(const CFreal resFactor);
  
  /// compute the Riemann flux to gradient jacobian numerically
  virtual void computeRiemannFluxToGradJacobianNum(const CFreal resFactor);
  
  /// compute the cell flux jacobian numerically
  virtual void computeCellFluxJacobianNum(const CFreal resFactor);
  
  /// compute the flux to gradient jacobian numerically
  virtual void computeFluxToGradJacobianNum(const CFreal resFactor);
  
  /**
   * compute the contribution of the diffusive face term to both Jacobians
   */
  virtual void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  virtual void computeOneJacobDiffFaceTerm(const CFuint side);
  
protected: // data

  /// handle to the wall distance
  Framework::DataSocketSink<CFreal> socket_wallDistance;
  
  /// idx of closest sol to each flx
  Common::SafePtr< std::vector< CFuint > > m_closestSolToFlxIdx;
  
  Common::SafePtr< Physics::KOmega::NavierStokes2DKLogOmega > m_navierStokesVarSetTurb;
  
  Common::SafePtr< Physics::KOmega::NavierStokes3DKLogOmega > m_navierStokesVarSetTurb3D;
  
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

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstructionTurb_hh
