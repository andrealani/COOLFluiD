// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/LLAVJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Command to add Localized Laplacian Artificial Viscosity near discontinuities 
 * for implicit schemes for NS
 * 
 * @author Ray Vandenhoeck
 */
class LLAVJacobFluxReconstructionNS : public LLAVJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit LLAVJacobFluxReconstructionNS(const std::string& name);

  /// Destructor
  virtual ~LLAVJacobFluxReconstructionNS() {}
  
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
   * Unsetup private data
   */
  virtual void unsetup();
  
protected: //functions

  /**
   * Set the data for the current face necessary to calculate FI
   */
  virtual void setFaceData(CFuint faceID);
  
  /**
   * Set the data for the current cell necessary to calculate the residual update
   */
  virtual void setCellData();
  
  /**
   * set the face neighbour gradients
   * @pre setOtherFacesLocalIdxs()
   * @pre setFaceNeighbourStates()
   */
  virtual void setFaceNeighbourGradients();


protected: //data
  
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionNS_hh

