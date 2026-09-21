// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_MFMHDBndGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_MFMHDBndGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"
#include "MultiFluidMHD/DiffMFMHDVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective boundary command for multi-fluid MHD with a diffusive term. The
 * boundary corrections of the gradients are computed by the base command.
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class MFMHDBndGradientComputer : public ConvBndCorrectionsRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit MFMHDBndGradientComputer(const std::string& name);

  /// Destructor
  virtual ~MFMHDBndGradientComputer() {}
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_MFMHDBndGradientComputer_hh

