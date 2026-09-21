// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_PoissonJacobBndGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_PoissonJacobBndGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstruction.hh"
#include "Poisson/PoissonDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective boundary command for implicit schemes for Poisson. The Poisson
 * equation has no convective flux, so this command only adds the boundary face
 * corrections of the gradients, computed by the base command.
 * 
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class PoissonJacobBndGradientComputer : public ConvBndCorrectionsRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit PoissonJacobBndGradientComputer(const std::string& name);

  /// Destructor
  virtual ~PoissonJacobBndGradientComputer() {}

protected: //functions
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);
  
  /**
   * Add the boundary face corrections of the gradients for the faces of the
   * current TRS. No boundary flux, wave speed or Jacobian contribution is
   * computed, since the Poisson equation has no convective flux.
   */
  virtual void executeOnTrs();
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PoissonJacobBndGradientComputer_hh

