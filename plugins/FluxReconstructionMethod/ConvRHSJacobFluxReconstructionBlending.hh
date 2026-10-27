// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstructionBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstructionBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/BlendingP0Data.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Convective RHS and Jacobian command with cell-wise P0 order
 * blending (UseBlending = true), implicit counterpart of
 * ConvRHSFluxReconstructionBlending.
 *
 * Notation: alpha_K the blending coefficient of cell K (socket "alpha"),
 * alpha_f = max(alpha_L, alpha_R) the coefficient of a face, F_I the FR interface
 * flux, F* the same Riemann flux evaluated with the conserved cell averages of
 * the two neighbours, W_fp the cell-average weight of the correction of flux
 * point fp (BlendingP0Data), divh_i,fp the divergence of the correction
 * function of fp at solution point i. Residual of cell K at solution point i:
 *   R_i = (1 - alpha_K) [ -div F_D + sum_fp divh_i,fp F_D - sum_fp divh_i,fp F_I ]
 *       - sum_fp W_fp [ alpha_f F*_fp + (alpha_K - alpha_f) F_I ]
 * The last sum is the same at every solution point and moves the cell average
 * by -sum_fp W_fp [(1 - alpha_f) F_I + alpha_f F*], a flux both neighbours
 * share for any alpha field.
 *
 * Jacobian, alpha frozen: the volume term is (1 - alpha_K) times the FR volume
 * Jacobian. On a face with alpha_f > 0 every solution point of both cells
 * depends on every solution point of both cells (cell averages in F* and the
 * term that is the same at every solution point), so the face block is dense;
 * faces with alpha_f = 0 use the FR face Jacobian.
 *
 * @author Rayan Dhib
 */
class ConvRHSJacobFluxReconstructionBlending : public ConvRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvRHSJacobFluxReconstructionBlending(const std::string& name);

  /// Destructor
  virtual ~ConvRHSJacobFluxReconstructionBlending() {}

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Execute processing actions: residual, then the face and cell Jacobians
  virtual void execute();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// FR interface flux, and the P0 flux of the cell averages when alpha_f > 0
  virtual void computeInterfaceFlxCorrection();

  /// FR face correction scaled by (1 - alpha_K) plus the mean correction of the face
  virtual void computeCorrection(CFuint side, std::vector< RealVector >& corrections);

  /// FR volume term scaled by (1 - alpha_K)
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);

  /// FR volume term of the perturbed state scaled by (1 - alpha_K)
  virtual void computePertDivDiscontFlx(RealVector& residuals);

  /// flux projection vectors and alpha of the current cell
  virtual void setCellData();

  /// FR wave speed scaled by (1 - alpha) plus the P0 wave speed scaled by alpha
  virtual void updateWaveSpeed();

  /// blending coefficient of the cell holding the given states
  CFreal getCellAlpha(const std::vector< Framework::State* >& states);

  /// cell averages of both neighbours of the current face
  void computeFaceCellAverages();

  /// P0 flux F* of the cell averages at every flux point of the current face
  void computeFaceFluxP0();

  /// P0 wave speed of one side of the current face, from the cell average of that side
  CFreal computeWaveSpeedP0(const CFuint side);

  /// mean correction -sum_fp W_fp [alpha_f F* + (alpha_K - alpha_f) F_I], added to every solution point of one side
  void addMeanFaceTerm(const CFuint side, const CFreal alpha, std::vector< RealVector >& corrections);

  /**
   * Face part of the residual of one side at every solution point, from the
   * current F_I (m_cellFlx), cell averages and F*, as a flat vector
   * [solution point][equation]. Used for the unperturbed and the perturbed
   * values of the dense face Jacobian, so both go through the same operations.
   */
  void computeFullFaceCorrection(const CFuint side, RealVector& corrections);

  /**
   * Dense face Jacobian for alpha_f > 0.
   * @param side  side whose rows are filled; both sides when bothSides is true
   */
  void computeFaceJacobBlend(const bool bothSides, const CFuint side);

protected: // data

  /// socket holding the per-cell blending coefficient
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// weights and cell averages
  BlendingP0Data m_p0Data;

  /// cell averages of the two neighbours of the current face, update variables
  std::vector< Framework::State* > m_cellAvgStates;

  /// backup of the cell average of the perturbed side
  RealVector m_cellAvgBackup;

  /// true when the cell averages of the current face are up to date
  bool m_faceAvgsDone;

  /// weight of the P0 flux at the current face
  CFreal m_currFaceAlphaF;

  /// blending coefficient of the current cell
  CFreal m_currCellAlpha;

  /// P0 flux F* of the cell averages in the flux points of the current face
  std::vector< RealVector > m_flxPntFluxP0;

  /// backup of m_flxPntFluxP0 during the perturbations
  std::vector< RealVector > m_flxPntFluxP0Backup;

  /// mean correction of one side of the current face, per equation
  RealVector m_meanTerm;

}; // class ConvRHSJacobFluxReconstructionBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstructionBlending_hh
