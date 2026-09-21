// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_MHDGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_MHDGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstruction.hh"
#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective command for MHD with a diffusive term. The gradients of the
 * gradient variables are computed by the base command.
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class MHDGradientComputer : public ConvRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit MHDGradientComputer(const std::string& name);

  /// Destructor
  virtual ~MHDGradientComputer() {}

}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_MHDGradientComputer_hh

