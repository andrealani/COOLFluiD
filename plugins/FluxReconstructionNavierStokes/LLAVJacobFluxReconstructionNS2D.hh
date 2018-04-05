// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionNS2D_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionNS2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/LLAVJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// Command to add Localized Laplacian Artificial Viscosity near discontinuities for implicit schemes
/// @author Ray Vandenhoeck
    
class LLAVJacobFluxReconstructionNS2D : public LLAVJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit LLAVJacobFluxReconstructionNS2D(const std::string& name);

  /// Destructor
  virtual ~LLAVJacobFluxReconstructionNS2D() {}
  
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
  
  /**
   * set the necessary variables to compute gradients
   */
  virtual void setNecGradientVars(const std::vector<RealVector*>& states, RealMatrix& values, const CFuint stateSize);

protected: //data
  
  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionNS2D_hh

