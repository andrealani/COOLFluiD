// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the system using a FluxReconstruction solver
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class DiffRHSFluxReconstructionNS : public DiffRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit DiffRHSFluxReconstructionNS(const std::string& name);

  /// Destructor
  virtual ~DiffRHSFluxReconstructionNS() {}

protected: //functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * Add the updateCoeff corrections due to the partition faces
   */
  void addPartitionFacesCorrection();
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSFluxReconstructionNS_hh

