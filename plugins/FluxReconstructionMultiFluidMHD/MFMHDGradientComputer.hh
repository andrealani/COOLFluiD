// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_MFMHDGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_MFMHDGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstruction.hh"
#include "MultiFluidMHD/DiffMFMHDVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective command for multi-fluid MHD with a diffusive term. The gradients
 * of the gradient variables are computed by the base command.
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class MFMHDGradientComputer : public ConvRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit MFMHDGradientComputer(const std::string& name);

  /// Destructor
  virtual ~MFMHDGradientComputer() {}
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffMFMHDVarSetGradientComputer_hh

