// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstructionMHD_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstructionMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvDiffJacobFluxReconstruction.hh"
#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes the contribution of the 
 * convective, diffusive to the RHS for a Flux Reconstruction scheme for implicit 
 * time marching for MHD
 * 
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class ConvDiffJacobFluxReconstructionMHD : public ConvDiffJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvDiffJacobFluxReconstructionMHD(const std::string& name);

  /// Destructor
  virtual ~ConvDiffJacobFluxReconstructionMHD() {}
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: //functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
protected: // data
  
  /// diffusive variable set
  Common::SafePtr< Physics::MHD::MHDProjectionDiffVarSet > m_diffusiveVarSetMHD;
  
  /// Vector transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> m_updateToSolutionVecTrans;
  
  /// number of species
  CFuint m_nbrSpecies;
  
  /// unperturbed grad vars
  RealVector m_unpertGradVars;
  
  /// perturbed grad vars
  RealVector m_pertGradVars;
  
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

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstructionMHD_hh
