// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes the contribution of the 
 * diffusive tems to the RHS for a Flux Recosntruction scheme for implicit 
 * time marching for NS
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
class DiffRHSJacobFluxReconstructionNS : public DiffRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit DiffRHSJacobFluxReconstructionNS(const std::string& name);

  /// Destructor
  virtual ~DiffRHSJacobFluxReconstructionNS() {}

protected: //functions

  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * compute the terms for the gradient computation for a bnd face
   */
  virtual void computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm);
  
  /**
   * compute the term for the gradient computation for the cell
   */
  virtual void computeCellGradTerm(RealMatrix& gradTerm);
  
  /**
   * compute the terms for the gradient computation for a face
   */
  virtual void computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR);
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
protected: // data
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionNS_hh

