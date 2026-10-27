// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/BlendingP0Data.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Explicit convective RHS command with cell-wise P0 order blending
 * (UseBlending = true).
 *
 * Notation: alpha_K the blending coefficient of cell K (socket "alpha", filled by
 * the OrderBlending command), alpha_f = max(alpha_L, alpha_R) the coefficient of a
 * face, F_D the discontinuous flux of the cell extrapolated to its flux points, F_I
 * the FR interface Riemann flux, F* the same Riemann flux evaluated with the
 * conserved cell averages u_bar_L and u_bar_R, divh_i,fp the divergence of the
 * correction function of flux point fp at solution point i, W_fp the cell-average
 * weight of that correction (BlendingP0Data).
 *
 * Residual of cell K at solution point i:
 *   R_i = (1 - alpha_K) [ -div F_D + sum_fp divh_i,fp F_D - sum_fp divh_i,fp F_I ]   (FR part)
 *       - sum_fp W_fp [ alpha_f F*_fp + (alpha_K - alpha_f) F_I,fp ]                 (low-order part)
 * The low-order part is the same at every solution point, so it changes the cell
 * average only. It equals alpha_K R_P0 - sum_fp W_fp (alpha_f - alpha_K)(F* - F_I),
 * with R_P0 = -sum_fp W_fp F* the first-order finite-volume residual of the cell
 * averages; the second term acts only on faces where the neighbour is more blended.
 * The cell average moves by -sum_fp W_fp [(1 - alpha_f) F_I + alpha_f F*], a flux
 * both neighbours share for any alpha field.
 * alpha = 0 everywhere is plain FR. Boundary faces use the FR boundary flux
 * (ConvBndCorrectionsRHSFluxReconstructionBlending).
 *
 * @author Rayan Dhib
 */
class ConvRHSFluxReconstructionBlending : public ConvRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvRHSFluxReconstructionBlending(const std::string& name);

  /// Destructor
  virtual ~ConvRHSFluxReconstructionBlending() {}

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// FR interface flux, and the P0 flux of the cell averages when alpha_f > 0
  virtual void computeInterfaceFlxCorrection();

  /// FR face correction scaled by (1 - alpha_K) plus the low-order part of the face
  virtual void computeCorrection(CFuint side, std::vector< RealVector >& corrections);

  /// FR volume term scaled by (1 - alpha)
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);

  /// flux projection vectors and alpha of the current cell
  virtual void setCellData();

  /// FR wave speed scaled by (1 - alpha) plus the P0 wave speed scaled by alpha
  virtual void updateWaveSpeed();

  /// blending coefficient of the cell holding the given states
  CFreal getCellAlpha(const std::vector< Framework::State* >& states);

  /// cell averages of both neighbours of the current face
  void computeFaceCellAverages();

  /// P0 wave speed of one side of the current face, from the cell average of that side
  CFreal computeWaveSpeedP0(const CFuint side);

  /// low-order part of the face, -sum_fp W_fp [alpha_f F* + (alpha_K - alpha_f) F_I],
  /// added to every solution point of one side
  void addMeanFaceTerm(const CFuint side, const CFreal alpha, std::vector< RealVector >& corrections);

protected: // data

  /// socket holding the per-cell blending coefficient
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// weights and cell averages
  BlendingP0Data m_p0Data;

  /// cell averages of the two neighbours of the current face, update variables
  std::vector< Framework::State* > m_cellAvgStates;

  /// true when the cell averages of the current face are up to date
  bool m_faceAvgsDone;

  /// weight of the P0 flux at the current face
  CFreal m_currFaceAlphaF;

  /// blending coefficient of the current cell
  CFreal m_currCellAlpha;

  /// P0 flux F* of the cell averages in the flux points of the current face
  std::vector< RealVector > m_flxPntFluxP0;

}; // class ConvRHSFluxReconstructionBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionBlending_hh
