// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstructionPoisson_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstructionPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/ConvDiffJacobFluxReconstructionNS.hh"

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
 * convective, diffusive to the RHS for a Flux Reconstruction scheme for implicit 
 * time marching for Poisson
 * 
 * @author Ray Vandenhoeck
 */
class ConvDiffJacobFluxReconstructionPoisson : public ConvDiffJacobFluxReconstructionNS {

public: // functions

  /// Constructor
  explicit ConvDiffJacobFluxReconstructionPoisson(const std::string& name);

  /// Destructor
  virtual ~ConvDiffJacobFluxReconstructionPoisson() {}
  
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
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected: //functions
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /// compute the Riemann flux to gradient jacobian numerically
  virtual void computeRiemannFluxToGradJacobianNum(const CFreal resFactor);
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * Compute the discontinuous contribution to the corrected gradients
   */
  virtual void computeGradients();
  
  /**
   * Compute the correction part of the corrected gradient
   */
  virtual void computeGradientFaceCorrections();
  
  /// compute the Riemann flux jacobian numerically
  virtual void computeRiemannFluxJacobianNum(const CFreal resFactor);
  
  /// compute the gradient variables to state jacobians numerically
  virtual void computeGradVarsToStateJacobianNum();
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
protected: // data
  
  Common::SafePtr<Physics::Poisson::PoissonDiffVarSet> m_diffVarSetPoisson;
  Common::SafePtr<Physics::Poisson::PoissonConvVarSet> m_convVarSetPoisson;
  
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

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffJacobFluxReconstructionPoisson_hh
