// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSGradientComputerSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_NSGradientComputerSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective command for NS with a diffusive term and subcell blending. The
 * gradients of the gradient variables are computed by the base command.
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class NSGradientComputerSubcellBlending : public ConvRHSFluxReconstructionSubcellBlending {

public: // functions

  /// Constructor
  explicit NSGradientComputerSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~NSGradientComputerSubcellBlending() {}
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSGradientComputerSubcellBlending_hh
