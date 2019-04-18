// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVDiffFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVDiffFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/LLAVDiffFluxReconstruction.hh"

#include "Framework/MultiScalarTerm.hh"

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
      class EulerTerm;
      class Euler2DVarSet;
    }
  }
  
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Command to add Localized Laplacian Artificial Viscosity near discontinuities 
 * for NS
 * 
 * @author Ray Vandenhoeck
 */
class LLAVDiffFluxReconstructionNS : public LLAVDiffFluxReconstruction {

public: // functions

  /// Constructor
  explicit LLAVDiffFluxReconstructionNS(const std::string& name);

  /// Destructor
  virtual ~LLAVDiffFluxReconstructionNS() {}
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options)
  {
  }
  
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
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * Compute the Peclet number based on the user input
   */
  virtual CFreal computePeclet();
  
  /**
   * command to compute the coefficient to be multiplied with epsilon for the wavespeedupdates
   */
  virtual CFreal computeViscCoef(RealVector* state);
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness();
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /**
   * Compute the left and right states and gradients in the flx pnts
   */
  virtual void computeFlxPntStatesAndGrads();
  
  /// compute the divergence of the discontinuous flux (-divFD+divhFD)
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);

protected: //data
  
  /// backup of the gradients in the neighbouring cell
  std::vector< std::vector< std::vector< RealVector* > > > m_gradsBackUp;
  
  /// physical model
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_eulerVarSet;
  
  /// physical model MS
  Common::SafePtr< Framework::MultiScalarTerm< Physics::NavierStokes::EulerTerm > > m_msEulerTerm;
  
  /// number of species
  CFuint m_nbrSpecies;
  
  /// physical model (in conservative variables)
  Common::SafePtr< Physics::NavierStokes::MultiScalarVarSet< Physics::NavierStokes::Euler2DVarSet > > m_eulerVarSet2;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T) for left neighbor
  RealMatrix m_tempGradTermL;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T) for right neighbor
  RealMatrix m_tempGradTermR;
  
  /// diffusive variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;
  
  /// element states of the left neighbor in the correct format
  std::vector< RealVector* > m_tempStatesL;
  
  /// element states of the right neighbor in the correct format
  std::vector< RealVector* > m_tempStatesR;
  
  /// the gradients in the neighbouring cell
  std::vector< std::vector< std::vector< RealVector >* > > m_cellGradsAV;
  
  /// average gradients in a flux point
  std::vector< RealVector* > m_avgGradAV;
  
  /// the corrected gradients in the flux points
  std::vector< std::vector< std::vector< RealVector* > > > m_cellGradFlxPntAV;
  
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

#endif // COOLFluiD_FluxReconstructionMethod_LLAVDiffFluxReconstructionNS_hh

