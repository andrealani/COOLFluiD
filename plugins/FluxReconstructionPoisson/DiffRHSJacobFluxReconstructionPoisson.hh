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
 * @author Rayan Dhib
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
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
  /**
   * Compute the unperturbed diffusive residual of the cell on one side of the
   * current face with the base command, then store the magnetic field of that
   * cell in the B sockets.
   * @param side side of the cell, LEFT or RIGHT
   */
  virtual void computeUnpertCellDiffResiduals(const CFuint side);

  /**
   * Store the magnetic field B = grad(phi) at the solution points of the cell on
   * one side of the current face, with phi the potential and grad(phi) its
   * corrected gradient, in the sockets Bx, By, Bz and, with (x,y,z) the
   * coordinates of the solution point, r = sqrt(x^2+y^2+z^2) and
   * rXY = sqrt(x^2+y^2), in the sockets Br = (x*Bx + y*By + z*Bz)/r,
   * Btheta = -y*Bx + x*By and Bphi = (z*x*Bx + z*y*By)/rXY - rXY*Bz.
   * @param side side of the cell, LEFT or RIGHT
   */
  void computeMagneticField(const CFuint side);
  
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

