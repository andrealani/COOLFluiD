// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionGammaAlpha_hh
#define COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionGammaAlpha_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionTurb/DiffBndCorrectionsRHSFluxReconstructionTurb.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes the contribution of the
   * boundary faces of the diffusive terms to the RHS for the Gamma-Alpha
   * transition model. It is the turbulence command
   * DiffBndCorrectionsRHSFluxReconstructionTurb under its own name: the
   * transition flag the Gamma-Alpha wall boundary condition reads is set once
   * per face in the gradient stage (GammaAlphaBndGradientComputer) and kept on
   * the boundary condition.
   *
   * @author Ray Vandenhoeck
   * @author Rayan Dhib
   */
class DiffBndCorrectionsRHSFluxReconstructionGammaAlpha : public DiffBndCorrectionsRHSFluxReconstructionTurb {

public:

  /**
   * Constructor
   */
  DiffBndCorrectionsRHSFluxReconstructionGammaAlpha(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffBndCorrectionsRHSFluxReconstructionGammaAlpha();

}; // end of class DiffBndCorrectionsRHSFluxReconstructionGammaAlpha

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_DiffBndCorrectionsRHSFluxReconstructionGammaAlpha_hh
