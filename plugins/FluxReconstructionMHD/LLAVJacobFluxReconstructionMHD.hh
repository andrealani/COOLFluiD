// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionMHD_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/LLAVJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }
  
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Command to add Localized Laplacian Artificial Viscosity near discontinuities 
 * for implicit schemes for MHD
 * 
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class LLAVJacobFluxReconstructionMHD : public LLAVJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit LLAVJacobFluxReconstructionMHD(const std::string& name);

  /// Destructor
  virtual ~LLAVJacobFluxReconstructionMHD() {}
  
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
   * Compute the Peclet number based on the user input
   */
  virtual CFreal computePeclet();
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness();
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness(const CFuint side);

protected: //data
  
  /// physical model
  Common::SafePtr< Physics::MHD::MHD3DProjectionVarSet > m_varSet;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  /// Physical data temporary vector
  RealVector m_pData2;
  
    /// Inverse of the Vandermonde
    RealMatrix m_vdmInv;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionMHD_hh

