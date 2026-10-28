// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstructionBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstructionBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/BlendingP0Data.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Boundary face corrections of the cell-wise P0 order blending
 * (UseBlending = true), companion of ConvRHSFluxReconstructionBlending.
 *
 * The face flux is the FR boundary flux Fhat_fp of the parent class (no boundary
 * condition is evaluated on cell averages, and periodic faces see the same flux
 * from both sides). With alpha the blending coefficient of the cell, W_fp the
 * cell-average weight of the correction of flux point fp (BlendingP0Data)
 * and divh_i,fp the correction function divergence, the correction of solution
 * point i is
 *   -(1 - alpha) sum_fp divh_i,fp Fhat_fp - alpha sum_fp W_fp Fhat_fp.
 *
 * @author Rayan Dhib
 */
class ConvBndCorrectionsRHSFluxReconstructionBlending : public ConvBndCorrectionsRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvBndCorrectionsRHSFluxReconstructionBlending(const std::string& name);

  /// Destructor
  virtual ~ConvBndCorrectionsRHSFluxReconstructionBlending() {}

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// FR boundary correction scaled by (1 - alpha) plus the P0 part of the face
  virtual void computeCorrection(std::vector< RealVector >& corrections);

  /// FR wave speed scaled by (1 - alpha) plus the P0 wave speed scaled by alpha
  virtual void updateWaveSpeed();

  /// blending coefficient of the current cell
  CFreal getCellAlpha();

  /// P0 wave speed of the current face, from the cell average
  CFreal computeWaveSpeedP0();

protected: // data

  /// socket holding the per-cell blending coefficient
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// metric Jacobian at the solution points, for the cell average
  Framework::DataSocketSink< CFreal > socket_volumes;

  /// weights and cell averages
  BlendingP0Data m_p0Data;

  /// cell average of the current cell, update variables
  Framework::State* m_cellAvgState;

  /// physical data of the cell average
  RealVector m_pDataP0;

}; // class ConvBndCorrectionsRHSFluxReconstructionBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstructionBlending_hh
