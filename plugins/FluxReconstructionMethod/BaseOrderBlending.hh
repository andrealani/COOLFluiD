// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_BaseOrderBlending_hh
#define COOLFluiD_FluxReconstructionMethod_BaseOrderBlending_hh

//////////////////////////////////////////////////////////////////////////////

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
 * Physics-agnostic base class. Handles monitored expressions that can be
 * extracted from any state/physical data: rho, p, rho*p, p/rho, rho/p,
 * velocity_magnitude. Physics-specific expressions (e.g. B2 for MHD) are
 * added by overriding extractMonitoredField() in a subclass.
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

  /// One Jacobi smoothing iteration: reads from m_sweepSnapshot, writes to socket_alpha.
  /// alpha_new[i] = applyAlphaLimits(max(snapshot[i], NeighborWeight * max_{j in N(i)} snapshot[j]))
  void applyJacobiSmoothingPass();

protected: // data

  /// Blending coefficient per solution point, consumed by the blending RHS classes.
  Framework::DataSocketSource< CFreal > socket_alpha;

  /// Previous-iteration alpha, used for the freeze mechanism.
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

  /// Minimum blending coefficient (dead-band at 0 and 1).
  CFreal m_alphaMin;

  /// Maximum blending coefficient cap.
  CFreal m_alphaMax;

  /// Decay factor applied to neighbor alpha during Jacobi smoothing.
  CFreal m_neighborWeight;

  /// Number of Jacobi smoothing iterations beyond the initial spread.
  /// Total spreading passes = m_nbSweeps + 1.
  CFuint m_nbSweeps;

  /// Iteration number at which alpha is frozen (reuses prevAlpha).
  CFuint m_freezeFilterIter;

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

}; // class BaseOrderBlending

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BaseOrderBlending_hh
