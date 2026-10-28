// Copyright (C) 2026 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_MHDJacobBlendingCommands_hh
#define COOLFluiD_FluxReconstructionMethod_MHDJacobBlendingCommands_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstructionBlending.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstructionBlending.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Order-blended convective command for implicit schemes for MHD with a
 * diffusive term: the blended counterpart of ConvRHSJacobMHD, selected by the
 * solver when UseBlending is on and ConvSolveCom = ConvRHSJacobMHD. The
 * gradients for the diffusive command are computed by the base command.
 *
 * @author Rayan Dhib
 */
class MHDJacobGradientComputerBlending : public ConvRHSJacobFluxReconstructionBlending {

public: // functions

  /// Constructor
  explicit MHDJacobGradientComputerBlending(const std::string& name);

  /// Destructor
  virtual ~MHDJacobGradientComputerBlending() {}

}; // class MHDJacobGradientComputerBlending

//////////////////////////////////////////////////////////////////////////////

/**
 * Order-blended convective boundary command for implicit schemes for MHD with
 * a diffusive term: the blended counterpart of ConvBndCorrectionsRHSJacobMHD,
 * selected by the solver when UseBlending is on and SpaceRHSJacobCom =
 * RHSJacobMHD.
 *
 * @author Rayan Dhib
 */
class MHDJacobBndGradientComputerBlending : public ConvBndCorrectionsRHSJacobFluxReconstructionBlending {

public: // functions

  /// Constructor
  explicit MHDJacobBndGradientComputerBlending(const std::string& name);

  /// Destructor
  virtual ~MHDJacobBndGradientComputerBlending() {}

}; // class MHDJacobBndGradientComputerBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_MHDJacobBlendingCommands_hh
