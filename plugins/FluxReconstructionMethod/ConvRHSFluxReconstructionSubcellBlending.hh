// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/SubcellBlendingQuadData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Explicit convective RHS command with subcell P0 order blending, for quads.
 *
 * At every solution point the residual is
 *   R = (1 - alpha) * R_FR + alpha * R_P0
 * where R_FR is the standard FR residual of the parent class and R_P0 is a
 * P0 residual on the subcell grid held by
 * SubcellBlendingQuadData. The per-cell blending coefficient alpha comes from
 * the "alpha" socket filled by the OrderBlending preprocessing command.
 *
 * Internal subcell interfaces use the configured Riemann flux between the two
 * adjacent solution point states.
 *
 * Element faces: the P0 boundary flux of the outer subcells is the face Riemann
 * flux already used by the FR correction, so both neighbours of a face exchange
 * the same flux and the blend is conservative for any alpha, also when alpha
 * differs between the two cells.
 *
 * FaceFluxBlending option (default true): the face Riemann flux is replaced by
 *   F* = (1 - alphaF) * F*(extrapolated states) + alphaF * F*(adjacent solution point states)
 * with alphaF = max(alpha_L, alpha_R). Both neighbours use the same alphaF, so
 * conservation is kept. With alpha = 1 the cell then runs a genuine P0
 * scheme on the subcells. Without this option the outer subcells of a cell
 * with Gauss-Legendre points are fed fluxes built from extrapolated states,
 * which are unlimited reconstructions.
 *
 * Companion commands: ConvBndCorrectionsRHSFluxReconstructionSubcellBlending for
 * boundary faces, and the Jacobian variants for implicit runs.
 *
 * @author Rayan Dhib
 */
class ConvRHSFluxReconstructionSubcellBlending : public ConvRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvRHSFluxReconstructionSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~ConvRHSFluxReconstructionSubcellBlending() {}

  /// Defines the Config Option's of this class
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures the command
  virtual void configure(Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected: // functions

  /// face Riemann flux, blended with the first-order flux between the adjacent
  /// solution points when FaceFluxBlending is on
  virtual void computeInterfaceFlxCorrection();

  /// FR face correction scaled by (1-alpha) plus the subcell P0 boundary flux scaled by alpha
  virtual void computeCorrection(CFuint side, std::vector< RealVector >& corrections);

  /// FR volume term scaled by (1-alpha) plus the internal subcell P0 fluxes scaled by alpha
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);

  /// metric terms of the current cell at the solution points and at the subcell interfaces
  virtual void setCellData();

  /// blending coefficient of the cell holding the given states
  CFreal getCellAlpha(const std::vector< Framework::State* >& states);

  /// weight of the first-order flux at the current face, max(alpha_L, alpha_R),
  /// zero when FaceFluxBlending is off
  CFreal computeFaceAlphaF();

  /// add the subcell P0 boundary flux of the current face to the corrections of one neighbour
  void addSubcellFaceFlux(const CFuint side, const CFreal alpha,
                          std::vector< RealVector >& corrections);

protected: // data

  /// socket holding the per-cell blending coefficient
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// subcell grid of the element
  SubcellBlendingQuadData m_scData;

  /// blend the face Riemann flux with the first-order flux from the adjacent solution points
  bool m_faceFluxBlending;

  /// weight of the first-order flux at the current face
  CFreal m_currFaceAlphaF;

  /// blending coefficient of the current cell
  CFreal m_currCellAlpha;

  /// local face index of each cell flux point
  Common::SafePtr< std::vector< CFuint > > m_flxPntFaceConn;

}; // class ConvRHSFluxReconstructionSubcellBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvRHSFluxReconstructionSubcellBlending_hh
