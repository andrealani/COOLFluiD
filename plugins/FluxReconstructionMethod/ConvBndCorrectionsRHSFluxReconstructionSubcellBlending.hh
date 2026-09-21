// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstructionSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstructionSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/SubcellBlendingQuadData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Boundary face contribution of the explicit convective RHS with subcell
 * P0 order blending, for quads. Companion of
 * ConvRHSFluxReconstructionSubcellBlending.
 *
 * The boundary Riemann flux, between the inner extrapolated state and the BC
 * ghost state, is used by both parts of the blend: the FR correction is scaled
 * by (1 - alpha) and the same flux is applied as the P0 boundary flux of the
 * subcell of the closest solution point, scaled by alpha.
 *
 * No low-order flux blending is done at boundary faces, because a low-order
 * variant would need the BC to produce a ghost state for a second set of inner
 * states.
 *
 * @author Rayan Dhib
 */
class ConvBndCorrectionsRHSFluxReconstructionSubcellBlending : public ConvBndCorrectionsRHSFluxReconstruction {

public: // functions

  /// Constructor
  ConvBndCorrectionsRHSFluxReconstructionSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~ConvBndCorrectionsRHSFluxReconstructionSubcellBlending();

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unset up private data and data of the aggregated classes in this command after processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// FR boundary correction scaled by (1-alpha) plus the subcell P0 boundary flux scaled by alpha
  virtual void computeCorrection(std::vector< RealVector >& corrections);

  /// blending coefficient of the cell behind the current boundary face
  CFreal getCellAlpha();

protected: // data

  /// socket holding the per-cell blending coefficient
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// subcell grid of the element
  SubcellBlendingQuadData m_scData;

}; // class ConvBndCorrectionsRHSFluxReconstructionSubcellBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstructionSubcellBlending_hh
