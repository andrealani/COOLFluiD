#ifndef COOLFluiD_FluxReconstructionMethod_BasePositivityIDP_hh
#define COOLFluiD_FluxReconstructionMethod_BasePositivityIDP_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/VarSetTransformer.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Physics-agnostic machinery for an invariant-domain-preserving (Zhang-Shu)
 * limiter for the FR solver.
 *
 * The limiter computes ONE scaling factor theta per cell, applied to the
 * conservative state about the metric-weighted cell mean, so that the
 * admissible-set constraints hold at every evaluation point. Because the
 * scaling is affine about the mean, the discretely conserved quantity
 * sum_j c_j |J_j| U_cons,j is preserved exactly.
 *
 * Evaluation point set = solution points + flux points + Zhang-Shu interior
 * nodes. The solution and flux points are what the scheme actually samples
 * (so they must be admissible or the wave speeds produce NaN); the Zhang-Shu
 * nodes are what the cell-mean update argument needs.
 *
 * NOTE on variables: the solver stores UPDATE variables, which need not be the
 * conservative set. When they are not, the scheme samples Cons(interp_Upd(x))
 * while this limiter scales the conservative nodal values, which defines
 * interp_Cons(x). Those differ, so the closed-form theta is a PREDICTOR and is
 * followed by a verify-and-bisect corrector; with conservative update variables
 * the two coincide and the corrector never fires. theta = 0 is always
 * admissible whatever the update set: it makes every nodal value equal, and the
 * Lagrange basis reproduces constants, so the interpolant is then constant and
 * equal to the cell mean.
 *
 * @author Rayan Dhib
 */

class BasePositivityIDP : public FluxReconstructionSolverCom {
public:

  /// which conservative components take part in the scaling
  enum ScaleMode { HYDRO = 0, FULL = 1, AUTO = 2 };

  explicit BasePositivityIDP(const std::string& name);

  virtual ~BasePositivityIDP();

  static void defineConfigOptions(Config::OptionList& options);

  virtual void setup();

  virtual void unsetup();

  virtual void configure ( Config::ConfigArgs& args );

  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /// Execute the limiter over all cells
  void execute();

protected: // physics hooks

  /**
   * Extract the constrained quantities from a conservative state.
   * @param cons  conservative state
   * @param rho   density
   * @param p     thermodynamic pressure
   * @param B2    squared magnetic field magnitude (0 for non-MHD physics)
   */
  virtual void constraintsAtPoint(const RealVector& cons,
                                  CFreal& rho, CFreal& p, CFreal& B2) const = 0;

  /// Map a conservative state back to update variables
  virtual void consToUpdate(const RealVector& cons, RealVector& update) const = 0;

  /// Conservative component indices scaled in the given mode
  virtual const std::vector< CFuint >& scaledIndices(ScaleMode mode) const = 0;

  /// Physics-specific setup, called at the end of setup(). Should validate
  /// the configured variable sets and throw if they are not supported.
  virtual void setupPhysics() = 0;

protected: // machinery

  /// Build the evaluation point set and the interpolation matrix (setup only)
  void buildEvalPointSet();

  /// Transform the stored update states of the current cell to conservative
  void computeCellConsStates();

  /// Metric-weighted conservative cell mean, using |J| from the volumes socket
  void computeMetricMean();

  /// Interpolate the conservative solution point values to the evaluation points
  void interpolateToEvalPnts(const std::vector< RealVector >& consSol,
                             std::vector< RealVector >& consEval) const;

  /// Interpolate the stored update states to the evaluation points, then
  /// transform to conservative. This is what the scheme actually samples.
  void computeTraceStates(std::vector< RealVector >& consEval);

  /// Build the Hydro-mode base state at a point into m_consBase: the cell mean
  /// for the scaled components, the local value for the frozen ones. This is
  /// the theta = 0 endpoint, so it is both the chord base for the pressure
  /// constraint and the test deciding the Auto fallback.
  void buildHydroBaseState(const RealVector& consPnt);

  /// theta for the density constraint at one evaluation point (exact, affine)
  CFreal thetaDensity(CFreal rhoBase, CFreal rho, CFreal epsRho) const;

  /// theta for the pressure constraint at one evaluation point (concavity chord)
  CFreal thetaPressure(CFreal pBase, CFreal p, CFreal epsP) const;

  /// Rebuild the conservative solution point values from the saved original
  /// state with the given scaling factor
  void applyThetaFromOrig(CFreal theta, ScaleMode mode);

  /// Scale the evaluation point states in place about the cell mean
  void scaleEvalStates(CFreal theta, ScaleMode mode);

  /// Write conservative solution point values back to the stored update states
  void writeBackStates();

  /// Check whether all evaluation points of consEval satisfy the constraints
  bool allAdmissible(const std::vector< RealVector >& consEval,
                     CFreal epsRho, CFreal epsP) const;

protected: // data

  /// socket for the positivity preservation feedback to LLAV
  Framework::DataSocketSink< CFreal > socket_posPrev;

  /// per solution point |J|, filled by StdSetup
  Framework::DataSocketSink< CFreal > socket_volumes;

  /// per solution point limiter activity, 1 - theta, for visualization
  Framework::DataSocketSource< CFreal > socket_outputPP;

  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  Framework::GeometricEntity* m_cell;

  std::vector< Framework::State* >* m_cellStates;

  /// update to solution (conservative) variable transformer
  Common::SafePtr< Framework::VarSetTransformer > m_updateToSolutionVecTrans;

  /// scratch state for the update to conservative transform
  Framework::State* m_tmpState;

  /// dummy coordinates for m_tmpState
  Framework::Node* m_tmpCoord;

  CFuint m_nbrEqs;

  CFuint m_dim;

  CFuint m_order;

  CFuint m_nbrSolPnts;

  /// total number of evaluation points (sol + flx + Zhang-Shu interior)
  CFuint m_nbrEvalPnts;

  /// number of Zhang-Shu interior nodes in the set, for reporting
  CFuint m_nbrZSPnts;

  /// number of points of the Gauss-Lobatto rule used by the Zhang-Shu argument
  CFuint m_nbrLobatto;

  /// first Gauss-Lobatto weight, normalized so the weights sum to one
  CFreal m_lobattoW1;

  /// interpolation from solution points to evaluation points
  std::vector< std::vector< CFreal > > m_evalPntPolyVals;

  /// normalized reference space quadrature weights
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;

  /// unlimited conservative states at the solution points of the current cell
  std::vector< RealVector > m_consSolOrig;

  /// working conservative states at the solution points of the current cell
  std::vector< RealVector > m_consSol;

  /// conservative states at the evaluation points of the current cell
  std::vector< RealVector > m_consEval;

  /// conservative cell mean
  RealVector m_consMean;

  /// scratch for the Hydro-mode base state
  RealVector m_consBase;

  /// scratch for a single update state
  RealVector m_tmpUpdate;

  /// configuration
  CFreal m_relDensityFloor;
  CFreal m_absDensityFloor;
  CFreal m_relPressureFloor;
  CFreal m_absPressureFloor;
  CFreal m_delta;
  CFuint m_maxVerifyIters;
  CFuint m_showrate;
  bool   m_enableDensity;
  bool   m_enablePressure;
  std::string m_scaleModeStr;
  ScaleMode m_scaleMode;

  /// diagnostics, per call
  CFuint m_nbLimited;
  CFuint m_nbMeanFailures;
  CFuint m_nbSequential;
  CFuint m_nbHydroFallbacks;
  CFuint m_nbVerifyBisections;
  CFuint m_nbVerifyToMean;
  CFreal m_minTheta;
  CFreal m_minRho;
  /// minimum pressure of the state arriving at this command, before limiting
  CFreal m_minP;
  /// minimum pressure of the state handed back to the scheme, after limiting.
  /// This one is bounded below by the floor by construction, so a negative
  /// value here is a bug, whereas a negative m_minP only says the time update
  /// produced an inadmissible state for the limiter to repair.
  CFreal m_minPOut;

}; // class BasePositivityIDP

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BasePositivityIDP_hh
