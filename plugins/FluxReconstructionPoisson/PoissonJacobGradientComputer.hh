// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_PoissonJacobGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_PoissonJacobGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstruction.hh"
#include "Poisson/PoissonDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective command for implicit schemes for Poisson. The Poisson equation has
 * no convective flux, so this command only computes the gradients of the
 * gradient variables with the functions of the base command.
 * 
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class PoissonJacobGradientComputer : public ConvRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit PoissonJacobGradientComputer(const std::string& name);

  /// Destructor
  virtual ~PoissonJacobGradientComputer() {}
  
  /**
   * Compute the gradients: the face corrections for every interior face and the
   * volume term for every cell. No interface flux, residual, wave speed or
   * Jacobian contribution is computed, since the Poisson equation has no
   * convective flux.
   */
  virtual void execute();
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PoissonJacobGradientComputer_hh

