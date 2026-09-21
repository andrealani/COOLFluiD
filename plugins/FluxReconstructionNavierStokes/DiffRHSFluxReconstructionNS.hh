// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffRHSFluxReconstruction.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes the contribution of the 
 * diffusive tems to the RHS for a Flux Recosntruction scheme for NS
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class DiffRHSFluxReconstructionNS : public DiffRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit DiffRHSFluxReconstructionNS(const std::string& name);

  /// Destructor
  virtual ~DiffRHSFluxReconstructionNS() {}
  
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
  
 protected: //data
  
  /// diffusive variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNS_hh

