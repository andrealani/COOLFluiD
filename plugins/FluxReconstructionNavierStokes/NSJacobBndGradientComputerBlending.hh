// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerBlending_hh
#define COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstructionBlending.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective boundary command for implicit schemes for NS with a diffusive
 * term when order blending is on. The boundary liftings of the gradients are
 * computed by the base command.
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class NSJacobBndGradientComputerBlending : public ConvBndCorrectionsRHSJacobFluxReconstructionBlending {

public: // functions

  /// Constructor
  explicit NSJacobBndGradientComputerBlending(const std::string& name);

  /// Destructor
  virtual ~NSJacobBndGradientComputerBlending() {}

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: //data
  
  /// scratch for derived commands: gradient variables at the flux points
  RealMatrix m_tempGradTerm;
  
  /// scratch for derived commands: gradient variables of the ghost states
  RealMatrix m_tempGradTermGhost;
  
  /// scratch for derived commands: flux point state data
  std::vector< RealVector* > m_tempStates;
  
  /// scratch for derived commands: ghost state data
  std::vector< RealVector* > m_tempStatesGhost;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerBlending_hh
