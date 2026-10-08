// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstructionSubcellBlending_hh
#define COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstructionSubcellBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/SubcellBlendingQuadData.hh"
#include "Common/SelfRegistPtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Implicit convective RHS command with subcell P0 order blending, for quads.
 *
 * Same blended residual as ConvRHSFluxReconstructionSubcellBlending, plus the
 * numerical Jacobian of that residual. One command serves both implicit modes:
 * with a full Jacobian assembly the perturbation loops run, in Jacobian-Free
 * Newton-Krylov mode the parent skips them and only the residual is used.
 *
 * The perturbed volume term uses a delta update. The unperturbed subcell
 * residual and the unperturbed interface fluxes are stored once per cell, and
 * perturbing a solution point only recomputes the interfaces that touch it, at
 * most two per reference direction, applying the difference. Both endpoints of
 * those interfaces share a reference line with the perturbed point, so nothing
 * is written outside the region the parent perturbed volume term fills.
 *
 * The sparsity is the same as without blending. The subcell boundary term at a
 * flux point touches only its closest solution point, which is already in the
 * flux point dependency list the face Jacobian accumulates over, and the
 * low-order face flux depends on that same solution point on each side.
 *
 * @author Rayan Dhib
 */
class ConvRHSJacobFluxReconstructionSubcellBlending : public ConvRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvRHSJacobFluxReconstructionSubcellBlending(const std::string& name);

  /// Destructor
  virtual ~ConvRHSJacobFluxReconstructionSubcellBlending();

  /// Defines the Config Option's of this class
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures the command
  virtual void configure(Config::ConfigArgs& args);

  /// records the reconstruction choices again when the blending hold stage changed (alpha
  /// frozen), then computes the residual (and the Jacobian) of the parent
  virtual void execute();

  /// Set up private data and data of the aggregated classes in this command before processing phase
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

  /// Returns the DataSocket's that this command provides as sources
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > providesSockets();

protected: // functions

  /// face Riemann flux, blended with the first-order flux between the adjacent
  /// solution points when FaceFluxBlending is on
  virtual void computeInterfaceFlxCorrection();

  /// FR face correction scaled by (1-alpha) plus the subcell P0 boundary flux scaled by alpha
  virtual void computeCorrection(CFuint side, std::vector< RealVector >& corrections);

  /// FR volume term scaled by (1-alpha) plus the internal subcell P0 fluxes scaled by alpha
  virtual void computeDivDiscontFlx(std::vector< RealVector >& residuals);

  /// true if the cell holding the given states uses the linear reconstruction:
  /// largest Mach number of its solution points below SubcellReconstructionMachMax
  /// record = true: unperturbed evaluation, the decision is stored while the blending is held
  bool reconstructCell(const std::vector< Framework::State* >& states, const CFuint cellID, const bool record);

  /// true once the order blending froze alpha (hold stage >= 1, see
  /// FluxReconstructionSolverData::getBlendingHoldStage)
  bool isReconstructionFrozen() const;

  /// hand the limiter freezing buffers to the subcell data, recording only if record and frozen
  void setLimiterRecording(const bool record);

  /// largest Mach number at the solution points of a cell, for SubcellReconstructionMachMax;
  /// physics-specific, 0 here (no restriction)
  virtual CFreal computeCellMaxMach(const std::vector< Framework::State* >& states) { return 0.; }

  /// store the element face samples of the current face for both cells (linear reconstruction)
  void storeFaceSamples();

  /// first-order flux at one element face flux point, between the closest solution point
  /// states or, with the linear reconstruction, between the reconstructed states
  const RealVector& computeFaceLoFlux(const CFuint iFlxPnt);

  /// metric terms of the current cell at the solution points and at the subcell interfaces
  virtual void setCellData();

  /// perturbed face Riemann flux of the influenced flux points, with the same blend
  virtual void computePertInterfaceFlxCorrection();

  /// perturbed FR face correction scaled by (1-alpha) plus the subcell P0 boundary flux
  virtual void computePertCorrection(CFuint side, RealVector& corrections);

  /// perturbed FR volume term scaled by (1-alpha) plus the subcell P0 delta update
  virtual void computePertDivDiscontFlx(RealVector& residuals);

  /// blending coefficient of the cell holding the given states
  CFreal getCellAlpha(const std::vector< Framework::State* >& states);

  /// weight of the first-order flux at the current face, max(alpha_L, alpha_R),
  /// zero when FaceFluxBlending is off
  CFreal computeFaceAlphaF();

protected: // data

  /// socket holding the per-cell blending coefficient
  Framework::DataSocketSink< CFreal > socket_alpha;

  /// element face samples for the linear reconstruction, [cell][flux point][equation]
  Framework::DataSocketSource< CFreal > socket_subcellFaceSamples;

  /// order of the subcell scheme at each solution point, for output (CGNS DataHandleOutput):
  /// 0 alpha = 0 (FR only), 1 first-order subcells, 2 linear reconstruction along both
  /// directions, 1.5 along one of them (the other set to first order by the admissibility test)
  Framework::DataSocketSource< CFreal > socket_subcellOrder;

  /// SubcellReconstruction option: FirstOrder or Linear
  std::string m_reconstruction;

  /// SubcellLimiter option: VanAlbada, Minmod or None
  std::string m_limiter;

  /// SubcellLimiterEps option: van Albada smoothing size, relative
  CFreal m_limiterEps;

  /// SubcellReconstructionVar option: variables of the reconstruction, empty for the update ones
  std::string m_reconstructionVar;

  /// SubcellReconstructionMachMax option: cells reaching this Mach number keep first-order subcells
  CFreal m_recMachMax;

  /// blending hold stage the stored choices belong to (0: nothing stored, live choices)
  CFuint m_recordedStage;

  /// true while the unperturbed face fluxes are computed (the mask may be recorded)
  bool m_recordFace;

  /// frozen reconstruction choice per cell: -1 not frozen yet, 0 first order, 1 reconstructed
  std::vector< CFint > m_cellMask;

  /// frozen limiter factors, [state][direction][equation], and their flags [state][direction]
  std::vector< CFreal > m_frozenPhi;
  std::vector< CFuint > m_frozenValid;

  /// transformers between the update and the reconstruction variables
  Common::SelfRegistPtr< Framework::VarSetTransformer > m_toRecTrans;
  Common::SelfRegistPtr< Framework::VarSetTransformer > m_fromRecTrans;

  /// reconstructed left and right states at an element face flux point
  std::vector< Framework::State* > m_faceRecStates;

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

  /// face flux of one flux point while it is being blended
  RealVector m_blendedFaceFlux;

}; // class ConvRHSJacobFluxReconstructionSubcellBlending

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvRHSJacobFluxReconstructionSubcellBlending_hh
