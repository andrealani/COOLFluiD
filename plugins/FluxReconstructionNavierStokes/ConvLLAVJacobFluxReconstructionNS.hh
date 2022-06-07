// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvLLAVJacobFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_ConvLLAVJacobFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvLLAVJacobFluxReconstruction.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

#include "Framework/MultiScalarTerm.hh"

#include "NavierStokes/MultiScalarVarSet.hh"

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
 * This class represents a command that computes the contribution of the 
 * convective, diffusive and LLAV tems to the RHS for a Flux Reconstruction scheme for implicit 
 * time marching for NS
 * 
 * @author Ray Vandenhoeck
 */
class ConvLLAVJacobFluxReconstructionNS : public ConvLLAVJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvLLAVJacobFluxReconstructionNS(const std::string& name);

  /// Destructor
  virtual ~ConvLLAVJacobFluxReconstructionNS() {}
  
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
  
  /// prepare the computation of the diffusive flux
  void prepareFluxComputation();
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /**
   * Compute the discontinuous contribution to the corrected gradients
   */
  virtual void computeGradients();
  
  /**
   * Compute the correction part of the corrected gradient
   */
  virtual void computeGradientFaceCorrections();
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness();
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness(const CFuint side);
  
  /**
   * Compute the reference artificial viscosity
   */
  virtual void computeEpsilon0();
  
  /**
   * Compute the reference artificial viscosity
   */
  virtual void computeEpsilon0(const CFuint side);
  
  /// compute the Riemann flux jacobian numerically
  virtual void computeRiemannFluxJacobianNum(const CFreal resFactor);
  
  /// compute the AV to state Jacobian
  virtual void computeEpsToStateJacobianAna();
  
protected: // data
  
  /// Vector transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> m_updateToSolutionVecTrans;
  
  /// physical model
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_eulerVarSet;
  
  /// physical model MS
  Common::SafePtr< Framework::MultiScalarTerm< Physics::NavierStokes::EulerTerm > > m_msEulerTerm;
  
  /// number of species
  CFuint m_nbrSpecies;
  
  /// physical model (in conservative variables)
  Common::SafePtr< Physics::NavierStokes::MultiScalarVarSet< Physics::NavierStokes::Euler2DVarSet > > m_eulerVarSet2;
  
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

#endif // COOLFluiD_FluxReconstructionMethod_ConvLLAVJacobFluxReconstructionNS_hh
