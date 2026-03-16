#ifndef COOLFluiD_FluxReconstructionMethod_EntropyFilter_hh
#define COOLFluiD_FluxReconstructionMethod_EntropyFilter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSource.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements an entropy-bounded exponential modal filter for the
 * Flux Reconstruction solver, following Dzanic & Witherden (JCP 2022).
 *
 * An exponential filter with factor f^(k^2) is applied per cell per mode
 * order k, where f = exp(-zeta) in [0,1]. The filter strength f is found
 * via bisection to enforce:
 *   (1) density positivity  rho > eps_rho
 *   (2) pressure positivity P   > eps_P
 *   (3) local minimum entropy principle  s = P * rho^{-gamma} >= s_min - e_tol
 *
 * Constraints are checked at both solution AND flux points.
 *
 * Physics-agnostic: uses ConvectiveVarSet::computePhysicalData() to extract
 * rho (index 0) and P (index 1) for any physics (Euler, NS, MHD).
 *
 * Hooked via the LimiterCom slot:  FluxReconstruction.LimiterCom = EntropyFilter
 *
 * @author Rayan Dhib
 */
class EntropyFilter : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit EntropyFilter(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~EntropyFilter();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Configures the command.
   */
  virtual void configure(Config::ConfigArgs& args);

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Execute Processing actions
   */
  void execute();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

protected: // helper functions

  /**
   * Compute filtered states at all solution points for a given filter parameter f.
   * Filtered state: u_tilde_j[iEq] = sum_{k=0}^{P} f^(k^2) * m_modeContrib[k][j][iEq]
   */
  void computeFilteredStates(const CFreal f);

  /**
   * Compute filtered state at a single evaluation point (sol or flx point).
   * Uses pre-decomposed mode contributions grouped by filter order.
   * @param ptIdx  index into m_modeContribPt arrays
   * @param f      filter parameter exp(-zeta)
   * @param filteredState  output filtered state
   */
  void computeFilteredStateAtPoint(const CFuint ptIdx, const CFreal f, RealVector& filteredState);

  /**
   * Compute density, pressure, and entropy at a state.
   * @param state  conservative state vector
   * @param rho    output density
   * @param p      output pressure
   * @param s      output specific entropy (p * rho^{-gamma})
   * @return true if state is physical (rho > 0 and p > 0)
   */
  bool computePhysics(const RealVector& state, CFreal& rho, CFreal& p, CFreal& s);

  /**
   * Pass 1: compute local minimum entropy per cell (sol + flx points).
   * Stores result in m_localMinEntropy[elemIdx].
   */
  void computeLocalMinEntropy();

  /**
   * Face loop: exchange entropy between face-adjacent cells.
   * For each internal face, both cells receive the min of both local entropies.
   * Result stored in m_neighborMinEntropy[elemIdx].
   */
  void exchangeNeighborEntropy();

protected: // data

  /// socket for output of the per-cell filter strength (1 - f_opt), 0 = unfiltered, 1 = max filter
  Framework::DataSocketSource< CFreal > socket_filterStrength;

  /// Heat capacity ratio (gamma)
  CFreal m_gamma;

  /// Minimum density floor
  CFreal m_minDensity;

  /// Minimum pressure floor
  CFreal m_minPressure;

  /// Maximum filter strength (zeta_max)
  CFreal m_zetaMax;

  /// Maximum root-bracketing iterations
  CFuint m_nBisect;

  /// Root-bracketing stopping tolerance
  CFreal m_tolerance;

  /// Entropy constraint tolerance (s >= sMin - eTol)
  CFreal m_eTol;

  /// Print frequency for filter statistics
  CFuint m_showrate;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of solution points for current element type
  CFuint m_nbrSolPnts;

  /// order of FR method
  CFuint m_order;

  /// Vandermonde matrix V (owned by frLocalData, persistent)
  Common::SafePtr< RealMatrix > m_vdm;

  /// Inverse Vandermonde V^{-1} (owned by frLocalData, persistent)
  Common::SafePtr< RealMatrix > m_vdmInv;

  /// Per-mode polynomial order (integer copy, avoids CFVec copy issues)
  std::vector< CFuint > m_maxModalOrder;

  /// Number of distinct mode orders (P+1 for tensor-product, may differ for simplices)
  CFuint m_nbrModeOrders;

  /// Coefficients for cell-averaged solution
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;

  /// Pre-decomposed mode contributions per order level
  /// m_modeContrib[k][iSol][iEq] = contribution of modes of order k to node iSol, eq iEq
  std::vector< std::vector< RealVector > > m_modeContrib;

  /// Filtered solution at sol points (working array)
  std::vector< RealVector > m_filteredStates;

  /// Cell-averaged state
  RealVector m_cellAvgState;

  /// Temporary state for computePhysicalData calls (allocated in setup, not constructor)
  Framework::State* m_tmpState;

  /// Dummy coordinate node for m_tmpState (some variable sets need coordinates)
  Framework::Node* m_tmpCoord;

  /// Physical data output vector
  RealVector m_physData;

  /// Access to the convective variable set (for computePhysicalData)
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// Modal coefficients working vector (nbrSolPnts)
  RealVector m_modalCoeffs;

  /// Modal coefficients in Legendre space (nbrSolPnts), working vector
  RealVector m_uHat;

  /// Number of flux points per cell
  CFuint m_nbrFlxPnts;

  /// Total number of evaluation points (sol + flx)
  CFuint m_nbrEvalPnts;

  /// Interpolation coefficients from sol points to flux points
  /// m_solPolyInFlxPnts[iFlxPnt][iSolPnt]
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyInFlxPnts;

  /// Pre-decomposed mode contributions at ALL evaluation points (sol + flx)
  /// m_modeContribPt[k][iPt] where iPt = 0..nbrSolPnts-1 for sol, nbrSolPnts..nbrEvalPnts-1 for flx
  std::vector< std::vector< RealVector > > m_modeContribPt;

  /// Filtered state at a single evaluation point (working array)
  RealVector m_filteredPtState;

  /// Face builder for internal face entropy exchange
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder>> m_faceBuilder;

  /// Per-cell local minimum entropy (own sol+flx points only)
  std::vector<CFreal> m_localMinEntropy;

  /// Per-cell neighbor-inclusive minimum entropy (after face exchange)
  std::vector<CFreal> m_neighborMinEntropy;

  /// Total number of cells across all element types
  CFuint m_totalCells;

  /// Number of cells filtered (local count)
  CFuint m_nbFiltered;

  /// Total number of cells filtered across all processors
  CFuint m_totalNbFiltered;

}; // class EntropyFilter

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_EntropyFilter_hh
