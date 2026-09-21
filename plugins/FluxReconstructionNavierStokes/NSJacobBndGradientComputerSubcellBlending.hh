// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective boundary command for implicit schemes for NS with a diffusive
 * term and subcell blending. The boundary liftings of the gradients are
 * computed by the base command.
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class NSJacobBndGradientComputerSubcellBlending : public ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending {

public: // functions

  /// Constructor
  explicit NSJacobBndGradientComputerSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~NSJacobBndGradientComputerSubcellBlending() {}
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSJacobBndGradientComputerSubcellBlending_hh
