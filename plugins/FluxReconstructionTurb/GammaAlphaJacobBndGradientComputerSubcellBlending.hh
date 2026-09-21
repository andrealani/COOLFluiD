// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_GammaAlphaJacobBndGradientComputerSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_GammaAlphaJacobBndGradientComputerSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/NSJacobBndGradientComputerSubcellBlending.hh"
#include "FluxReconstructionTurb/GammaAlphaTransitionCriterion.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective boundary command for implicit schemes with subcell order blending for
 * the Gamma-Alpha transition model: NSJacobBndGradientComputerSubcellBlending with the
 * transition flag of the wall boundary condition set before the ghost states
 * (GammaAlphaTransitionCriterion), once per face per residual evaluation.
 *
 * @author Rayan Dhib
 */
class GammaAlphaJacobBndGradientComputerSubcellBlending : public NSJacobBndGradientComputerSubcellBlending {

public: // functions

  /// Constructor
  explicit GammaAlphaJacobBndGradientComputerSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~GammaAlphaJacobBndGradientComputerSubcellBlending() {}

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: // functions

  /**
   * Sets the transition flags of the current face on the boundary condition.
   */
  virtual void prepareGhostStates();

protected: // data

  /// Jacobian determinant of every state
  Framework::DataSocketSink< CFreal > socket_volumes;

  /// metric of every state
  Framework::DataSocketSink< CFreal > socket_solPntNormals;

  /// the transition flags of a face
  GammaAlphaTransitionCriterion m_transitionCriterion;

}; // class GammaAlphaJacobBndGradientComputerSubcellBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_GammaAlphaJacobBndGradientComputerSubcellBlending_hh
