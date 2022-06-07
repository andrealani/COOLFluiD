// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionPoisson_hh
#define COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
    namespace Physics {
    namespace Poisson {
      class PoissonDiffVarSet;
      class PoissonConvVarSet;
    }
  }
    
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that computes the contribution of the 
 * diffusive tems to the RHS for a Flux Recosntruction scheme for implicit 
 * time marching for Poisson with numerical Jacobian
 * 
 * @author Ray Vandenhoeck
 */
class DiffRHSJacobFluxReconstructionPoisson : public DiffRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit DiffRHSJacobFluxReconstructionPoisson(const std::string& name);

  /// Destructor
  virtual ~DiffRHSJacobFluxReconstructionPoisson() {}
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

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
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /**
   * compute the unperturbed cell diffusive residuals
   * @pre m_faceTermComputers->computeDiffFaceTermAndUpdateCoefContributions
   * @pre setCellsData()
   */
  virtual void computeUnpertCellDiffResiduals(const CFuint side);
  
protected: // data
    
  /// matrix to store the state terms needed for the gradients (p, u, v, T) for left neighbor
  RealMatrix m_tempGradTermL;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T) for right neighbor
  RealMatrix m_tempGradTermR;
  
  Common::SafePtr<Physics::Poisson::PoissonDiffVarSet> m_diffVarSetPoisson;
  Common::SafePtr<Physics::Poisson::PoissonConvVarSet> m_convVarSetPoisson;
    
  /// element states of the left neighbor in the correct format
  std::vector< RealVector* > m_tempStatesL;
  
  /// element states of the right neighbor in the correct format
  std::vector< RealVector* > m_tempStatesR;
  
  /// damping coefficient
  CFreal m_dampCoeff;
  
  /// storage for Bx
  Framework::DataSocketSource<CFreal> socket_Bx;
  
  /// storage for By
  Framework::DataSocketSource<CFreal> socket_By;
  
  /// storage for Bz
  Framework::DataSocketSource<CFreal> socket_Bz;
  
  /// storage for Br
  Framework::DataSocketSource<CFreal> socket_Br;
  
  /// storage for Btheta
  Framework::DataSocketSource<CFreal> socket_Btheta;
  
  /// storage for Bphi
  Framework::DataSocketSource<CFreal> socket_Bphi;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffRHSJacobFluxReconstructionPoisson_hh

