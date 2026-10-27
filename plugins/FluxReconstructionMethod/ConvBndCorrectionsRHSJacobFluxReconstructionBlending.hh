// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/BlendingP0Data.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Boundary convective RHS and Jacobian command of the cell-wise P0
 * order blending (UseBlending = true), implicit counterpart of
 * ConvBndCorrectionsRHSFluxReconstructionBlending.
 *
 * The face flux at a boundary face is the FR boundary flux F_I (no P0 flux of
 * cell averages, alpha_f = 0). Residual of the cell K at solution point i, with
 * alpha_K its blending coefficient, divh_i,fp the divergence of the correction
 * function of flux point fp and W_fp its cell-average weight:
 *   C_i = -(1 - alpha_K) sum_fp divh_i,fp F_I,fp - alpha_K sum_fp W_fp F_I,fp
 * With alpha_K > 0 the second sum couples every solution point of the cell to
 * every one, so the Jacobian block of the cell is dense. alpha frozen.
 *
 * @author Rayan Dhib
 */
class ConvBndCorrectionsRHSJacobFluxReconstructionBlending : public ConvBndCorrectionsRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvBndCorrectionsRHSJacobFluxReconstructionBlending(const std::string& name);

  /// Destructor
  virtual ~ConvBndCorrectionsRHSJacobFluxReconstructionBlending() {}

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// FR boundary correction scaled by (1 - alpha_K) plus the P0 part -alpha_K sum_fp W_fp F_I
  virtual void computeCorrection(std::vector< RealVector >& corrections);

  /// perturbed correction: the full correction when alpha_K > 0, the FR one otherwise
  virtual void computePertCorrection(std::vector< RealVector >& corrections);

  /// every solution point is coupled to every one when alpha_K > 0
  virtual bool pertCouplesAllSolPnts()
  {
    return getCellAlpha() > 0.0;
  }

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

}; // class ConvBndCorrectionsRHSJacobFluxReconstructionBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSJacobFluxReconstructionBlending_hh
