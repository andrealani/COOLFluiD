// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/SubcellBlendingQuadData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Boundary face contribution of the implicit convective RHS with subcell
 * P0 order blending, for quads. Companion of
 * ConvRHSJacobFluxReconstructionSubcellBlending.
 *
 * The boundary Riemann flux, between the inner extrapolated state and the BC
 * ghost state, feeds both parts of the blend: the FR correction scaled by
 * (1 - alpha) and the P0 boundary flux of the subcell of the closest solution
 * point scaled by alpha. The same holds in the perturbed pass, so the numerical
 * Jacobian sees the blended residual.
 *
 * No low-order flux blending is done at boundary faces, because a low-order
 * variant would need the BC to produce a ghost state for a second set of inner
 * states.
 *
 * @author Rayan Dhib
 */
class ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending : public ConvBndCorrectionsRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending();

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unset up private data and data of the aggregated classes in this command after processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// FR boundary correction scaled by (1-alpha) plus the subcell P0 boundary flux scaled by alpha
  virtual void computeCorrection(std::vector< RealVector >& corrections);

  /// same blend, on the solution points affected by the perturbed one
  virtual void computePertCorrection(std::vector< RealVector >& corrections);

  /// blending coefficient of the cell behind the current boundary face
  CFreal getCellAlpha();

  /// add the subcell P0 boundary flux of one face flux point to the corrections
  void addSubcellFaceFlux(const CFuint iFlxPnt, const CFreal alpha, const CFint faceDir,
                          std::vector< RealVector >& corrections);

protected: // data

  /// socket holding the per-cell blending coefficient
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// subcell grid of the element
  SubcellBlendingQuadData m_scData;

  /// flags making sure each affected solution point is scaled exactly once
  std::vector< bool > m_scSolScaled;

}; // class ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending_hh
