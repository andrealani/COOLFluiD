// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputerSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputerSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective command for implicit schemes for NS with a diffusive term and
 * subcell blending. The gradients of the gradient variables are computed by
 * the base command.
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class NSJacobGradientComputerSubcellBlending : public ConvRHSJacobFluxReconstructionSubcellBlending {

public: // functions

  /// Constructor
  explicit NSJacobGradientComputerSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~NSJacobGradientComputerSubcellBlending() {}
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputerSubcellBlending_hh
