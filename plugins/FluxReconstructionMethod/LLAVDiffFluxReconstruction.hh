// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVDiffFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVDiffFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/LLAVFluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// Command to add Localized Laplacian Artificial Viscosity near discontinuities
/// @author Ray Vandenhoeck
/// @author Rayan Dhib
    
class LLAVDiffFluxReconstruction : public LLAVFluxReconstruction {

public: // functions

  /// Constructor
  explicit LLAVDiffFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~LLAVDiffFluxReconstruction() {}

  /// Execute processing actions
  void execute();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
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
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);

  /**
   * Compute the common flux at the flux points of the current face: the
   * physical diffusive flux with the compact gradients of the physical gradient
   * variables plus the artificial viscosity flux with the compact gradients of
   * the artificial viscosity variables.
   */
  virtual void computeInterfaceFlxCorrection();
  
  /**
   * Compute the divergence of the discontinuous flux (-divFD+divhFD) of the
   * current cell: the physical volume term of DiffRHSFluxReconstruction and
   * the artificial viscosity volume term and boundary flux of
   * LLAVFluxReconstruction, with the gradients switched to the artificial
   * viscosity gradients in between.
   */
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);

protected: //data

  /// physical common flux at the flux points of the current face
  std::vector< RealVector > m_physFlxPntRiemannFlux;

  /// physical divergence of the discontinuous flux at the solution points
  std::vector< RealVector > m_physDivContFlx;

  /// backup of the pointers to the gradients of the current cell
  std::vector< std::vector< RealVector >* > m_cellGradsPtrsBackUp;

  private:

  /// Physical data temporary vector
  RealVector m_pData;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LLAVDiffFluxReconstruction_hh

