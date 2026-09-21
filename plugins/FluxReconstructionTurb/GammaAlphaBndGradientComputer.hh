// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_GammaAlphaBndGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_GammaAlphaBndGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/NSBndGradientComputer.hh"
#include "FluxReconstructionTurb/GammaAlphaTransitionCriterion.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective boundary command for the Gamma-Alpha transition model:
 * NSBndGradientComputer with the transition flag of the wall boundary condition
 * set before the ghost states (GammaAlphaTransitionCriterion), once per face
 * per residual evaluation.
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class GammaAlphaBndGradientComputer : public NSBndGradientComputer {

public: // functions

  /// Constructor
  explicit GammaAlphaBndGradientComputer(const std::string& name);

  /// Destructor
  virtual ~GammaAlphaBndGradientComputer() {}

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

}; // class GammaAlphaBndGradientComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_GammaAlphaBndGradientComputer_hh
