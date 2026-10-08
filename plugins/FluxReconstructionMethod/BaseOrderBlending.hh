// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_BaseOrderBlending_hh
#define COOLFluiD_FluxReconstructionMethod_BaseOrderBlending_hh

//////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <deque>
#include <map>
#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Command that computes a per-cell order-blending coefficient alpha in [0,1]
 * based on a modal smoothness indicator (Persson-Peraire style), with
 * neighbor-max spreading via Jacobi smoothing passes.
 *
 * Produces sockets: alpha, prevAlpha, smoothness.
 *
 * Alpha is computed for every local cell, including the non-updatable overlap
 * cells of a parallel run, so that face-based blending reads a valid alpha
 * on both sides of a partition boundary.
 *
 * Physics-agnostic base class. Handles monitored expressions that can be
 * extracted from any state/physical data: rho, p, rho*p, p/rho, rho/p,
 * velocity_magnitude. Physics-specific expressions (e.g. B2 for MHD) are
 * added by overriding extractMonitoredField() in a subclass.
 *
 * Alpha can be frozen instead of following the sensor: from FreezeAlphaIter
 * on, or automatically once the residual stalls (FreezeAlphaStallWindow). The
 * subcell blending RHS freezes its reconstruction choices together with alpha
 * (FluxReconstructionSolverData::getBlendingHoldStage).
 *
 * @author Rayan Dhib
 */
class BaseOrderBlending : public FluxReconstructionSolverCom {
public:

  explicit BaseOrderBlending(const std::string& name);

  virtual ~BaseOrderBlending();

  virtual void setup();

  virtual void unsetup();

  virtual void configure(Config::ConfigArgs& args);

  static void defineConfigOptions(Config::OptionList& options);

  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  void execute();

  /// Vector of maximum modal order per mode for a given element shape and polynomial order.
  static RealVector getmaxModalOrder(const CFGeoShape::Type elemShape, const CFuint m_order);

protected: // functions

  /// Iteration from which alpha is frozen: the smaller of FreezeAlphaIter and its former
  /// name freezeFilterIter
  CFuint getFreezeAlphaIter() const { return std::min(m_freezeAlphaIter, m_freezeFilterIter); }

  /**
   * Fill m_tempSolPntVec with the monitored scalar evaluated at each
   * solution point of the current cell (m_cellStates).
   *
   * Base implementation handles: rho, p, rho*p, p/rho, rho/p, velocity_magnitude.
   * Physics-specific expressions must be added by overriding this method
   * in a physics-aware subclass (e.g. BaseOrderBlendingMHD for B2).
   */
  virtual void extractMonitoredField();

  /// Compute the per-cell modal smoothness indicator m_s (log10 of top-mode energy ratio).
  void computeSmoothness();

  /// Map log-smoothness to alpha via sinusoidal ramp in [s0-kappa, s0+kappa].
  /// s0 = -S0 * log10(order+1)
  CFreal computeBlendingCoefficient(CFreal smoothness) const;

  /// Apply the AlphaMin dead-band and AlphaMax cap to a raw alpha value.
  CFreal applyAlphaLimits(CFreal alpha) const;

  /**
   * Flag the unflagged cells that fail the undershoot test (hasUndershoot with
   * ForceAlphaMinMargin); flags are kept for the whole run. Releasing a flag
   * on a converged solution changes the operator at the front that needed it
   * (a flagged cell is first order, so the test cannot see whether it still
   * needs the flag). Does nothing when ForceAlphaMinVars is empty. Only owned
   * cells are tested (an overlap cell at the edge of the halo lacks some
   * neighbours); new flags are then sent to every rank, which updates its
   * copies of those cells.
   */
  void updateForcedCells();

  /// cell means of the ForceAlphaMinVars, and with ForceAlphaMinBoundaryGhost the smallest
  /// boundary face mean of each boundary cell, for all local cells (overlap included)
  void computeUndershootReferences();

  /// local face index of every boundary face (physical boundary TRSs, partition faces
  /// excluded) of each cell, built once from the boundary face lists of the solver data
  void buildBoundaryFaces();

  /// sets m_forcedCells to value on every rank's copy of the cells whose first state has one of
  /// the global IDs; returns the number of cells over all ranks
  CFuint shareForcedCells(const std::vector<CFuint>& firstStateGlobalIDs, const bool value);

  /// true if a ForceAlphaMinVars variable at a solution point of the current cell (m_cellStates)
  /// lies more than margin below the smallest cell mean of its neighbours and, with
  /// ForceAlphaMinBoundaryGhost, below the face means of its boundary faces
  bool hasUndershoot(const CFuint elemIdx, const CFreal margin);

  /**
   * Automatic freeze test (FreezeAlphaStallWindow = W > 0), called once per
   * iteration while alpha follows the sensor. Reads the monitored residual r of
   * the previous iteration and the CFL, both global, so every rank takes the
   * same decision: true once the CFL has not grown over the last W iterations
   * and r dropped by less than a tenth of a decade over them. Re-freezing later
   * from a more converged solution was tried and cycles: Mach-mask choices of
   * cells near SubcellReconstructionMachMax flip at every rebuild (Gnoffo P3).
   */
  bool isStalled();

  /// One Jacobi smoothing iteration: reads from m_sweepSnapshot, writes to socket_alpha.
  /// alpha_new[i] = applyAlphaLimits(max(snapshot[i], NeighborWeight * max_{j in N(i)} snapshot[j]))
  void applyJacobiSmoothingPass();

protected: // data

  /// Blending coefficient per solution point, consumed by the blending RHS classes.
  Framework::DataSocketSource< CFreal > socket_alpha;

  /// Previously applied alpha, used for temporal relaxation and freezing.
  Framework::DataSocketSource< CFreal > socket_prevAlpha;

  /// Per-solution-point smoothness indicator (for visualization / diagnostics).
  Framework::DataSocketSource< CFreal > socket_smoothness;

  /// Cell builder from the solver data.
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// Current cell geometric entity.
  Framework::GeometricEntity* m_cell;

  /// Pointer to the solution states of the current cell.
  std::vector< Framework::State* >* m_cellStates;

  /// Scratch buffer for Jacobi smoothing: pre-sweep snapshot of socket_alpha.
  std::vector< CFreal > m_sweepSnapshot;

  /// Update variable set, used by extractMonitoredField to access physical data.
  Common::SafePtr<Framework::ConvectiveVarSet> m_obUpdateVarSet;

  /// Physical data vector for extracting derived quantities (e.g. pressure).
  RealVector m_obPData;

  /// Expression selecting the monitored scalar field (e.g. "rho*p", "B2").
  std::string m_modalMonitoredExpression;

  /// Current element's temporary scalar values at solution points.
  RealVector m_tempSolPntVec;

  /// Scratch vector for the Vandermonde-inverse transform.
  RealVector m_tempSolPntVec2;

  /// Inverse Vandermonde matrix (nodal-to-modal transform).
  RealMatrix m_vdmInv;

  /// Per-mode maximum directional polynomial order (for energy ratio computation).
  RealVector m_maxModalOrder;

  /// Node-sharing neighbor IDs for each cell (pre-computed at setup).
  std::vector< std::vector<CFuint> > m_NeighborIDs;

  /// Current cell's smoothness indicator value (log10 of energy ratio).
  CFreal m_s;

  /// Reference smoothness threshold: s0 = -m_s0 * log10(order + 1).
  CFreal m_s0;

  /// Transition half-width around s0 for the sinusoidal ramp.
  CFreal m_kappa;

  /// Lower dead-band threshold: coefficients below this value become zero.
  CFreal m_alphaMin;

  /// Maximum blending coefficient cap.
  CFreal m_alphaMax;

  /// Decay factor applied to neighbor alpha during Jacobi smoothing.
  CFreal m_neighborWeight;

  /// Number of Jacobi smoothing iterations beyond the initial spread.
  /// Total spreading passes = m_nbSweeps + 1.
  CFuint m_nbSweeps;

  /// Iteration number at which alpha is frozen (reuses prevAlpha), option FreezeAlphaIter.
  CFuint m_freezeAlphaIter;

  /// Same, under the former option name freezeFilterIter; the smaller of the two is used.
  CFuint m_freezeFilterIter;

  /// Fraction of the newly computed sensor field applied after initialization; default 1 (no temporal relaxation).
  CFreal m_alphaRelaxation;

  /// True once prevAlpha holds an applied field for this setup/restart.
  bool m_alphaInitialized;

  /// FreezeAlphaStallWindow option: window W in iterations of the automatic freeze (0: off)
  CFuint m_freezeStallWindow;

  /// 0 while alpha follows the sensor, 1 once it is frozen (FreezeAlphaIter or FreezeAlphaStallWindow)
  CFuint m_stage;

  /// iteration of the last call, to test for a stall once per iteration
  CFuint m_lastIter;

  /// monitored residual and CFL of the last W + 1 iterations
  std::deque< CFreal > m_resHistory;
  std::deque< CFreal > m_cflHistory;

  /// FR polynomial order.
  CFuint m_order;

  /// Number of solution points per cell.
  CFuint m_nbrSolPnts;

  /// Number of equations in the physical model.
  CFuint m_nbrEqs;

  /// Dimensionality of the physical model.
  CFuint m_dim;

  /// Current element type index.
  CFuint m_iElemType;

  /// Current element index.
  CFuint m_elemIdx;

  /// state variables whose undershoot below the neighbouring cell means forces alpha = 1 (e.g. ln rho_i)
  std::vector< CFuint > m_forceAlphaMinVars;

  /// margin of the undershoot test on the ForceAlphaMinVars
  CFreal m_forceAlphaMinMargin;

  /// ForceAlphaMinBoundaryGhost option: boundary faces count as neighbours in the undershoot test
  bool m_forceAlphaMinBoundaryGhost;

  /// local face indices of the boundary faces of each cell (see buildBoundaryFaces)
  std::vector< std::vector< CFuint > > m_bndFaceOrients;

  /// true once m_bndFaceOrients is built
  bool m_bndFacesBuilt;

  /// flux points of each local face, and the solution polynomials at the flux points
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;

  /// cell means of the ForceAlphaMinVars, [cell][variable]
  std::vector< std::vector< CFreal > > m_minVarsCellMeans;

  /// smallest boundary face mean of the ForceAlphaMinVars, [cell][variable] (boundary cells)
  std::vector< std::vector< CFreal > > m_minVarsBndMeans;

  /// cells flagged by the undershoot test, alpha = 1 for the rest of the run
  std::vector< bool > m_forcedCells;

  /// local cell index of each cell, by the global ID of its first state (to apply the flags of other ranks)
  std::map< CFuint, CFuint > m_cellByFirstStateGlobalID;

}; // class BaseOrderBlending

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BaseOrderBlending_hh
