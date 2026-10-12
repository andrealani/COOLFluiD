#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_ManufacturedForcingTNEQ_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_ManufacturedForcingTNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNEQ/TNEQSourceTerm.hh"
#include "FluxReconstructionNEQ/Freno2021ManufacturedSolution.hh"

#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Forcing term of the manufactured solution of Freno et al. (2021), Eq. (19) (see Freno2021ManufacturedSolution),
 * for the inviscid two-temperature equations with the TNEQ source term. The steady discrete problem
 * -div F(U) + S(U) + f(x) = 0 has the manufactured state U_MS as exact solution when
 *   f(x) = div F(U_MS(x)) - S(U_MS(x)),
 * F being the convective flux and S the TNEQ source (chemistry and vibration-translation exchange). f is added at
 * the solution points like any source term. It is evaluated once per solution point, at the first residual:
 * - U_MS from the manufactured primitive state (RhoivtTv input variables) through the transformer to the update
 *   variables;
 * - div F by eighth-order central differences of the x and y fluxes of the update variable set, step h (FDStep), at
 *   x + k h e_x and x + k h e_y, k = -4..4, e_x and e_y the unit vectors;
 * - S by the TNEQSourceTerm code itself on temporary states holding U_MS.
 * f does not depend on the solution, so it has no Jacobian. Run the TNEQSourceTerm as a separate source term.
 * At the first residual the command logs, per equation, max |f|, max |S(U_MS)| and the FD check
 * max |divF_h - divF_2h| / max |divF_h|. Needs a dimensional 2D run (reference values 1, reference length 1).
 *
 * @author Rayan Dhib
 */
class ManufacturedForcingTNEQ : public TNEQSourceTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit ManufacturedForcingTNEQ(const std::string& name);

  /**
   * Destructor.
   */
  ~ManufacturedForcingTNEQ();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

protected:

  /**
   * Execute Processing actions: add the forcing, then log the checks after the first pass
   */
  void execute();

  /**
   * add the forcing at the solution points of the current cell
   */
  void addSourceTerm(RealVector& resUpdates);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// the forcing does not depend on the solution: zero Jacobian
  virtual void getSToStateJacobian(const CFuint iState);

  /// compute and store the forcing at the solution points of the current cell
  void computeForcingInCell();

  /**
   * manufactured state in update variables at a point
   * @param x, y coordinates of the point
   * @param state output state
   */
  void setManufacturedState(const CFreal x, const CFreal y, Framework::State& state);

  /**
   * divergence of the convective flux of the manufactured state by eighth-order central differences
   * @param x, y coordinates of the point
   * @param h finite difference step
   * @param divF output, one value per equation
   */
  void computeDivFlux(const CFreal x, const CFreal y, const CFreal h, RealVector& divF);

protected: // data

  /// option: species names in the order of the chemistry library
  std::vector<std::string> m_speciesNames;

  /// option: bar partial densities [kg/m3], same order as the names
  std::vector<CFreal> m_rhoBar;

  /// option: velocity scale ub [m/s]
  CFreal m_uBar;

  /// option: temperature scale Tb [K]
  CFreal m_TBar;

  /// option: vibrational temperature scale Tvb [K]
  CFreal m_TvBar;

  /// option: amplitude e of the perturbations
  CFreal m_eps;

  /// option: finite difference step h [m]
  CFreal m_fdStep;

  /// manufactured solution
  Freno2021ManufacturedSolution m_mms;

  /// transformer from the RhoivtTv input variables to the update variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;

  /// manufactured primitive state (rho_s, u, v, T, Tv)
  RealVector m_prim;

  /// input state for the transformer
  Framework::State* m_inputState;

  /// state used for the flux evaluations of the finite differences
  Framework::State* m_fdState;

  /// manufactured states at the solution points of a cell
  std::vector<Framework::State*> m_msStates;

  /// list of pointers handed to the TNEQ source term in place of the cell states
  std::vector<Framework::State*> m_msCellStates;

  /// physical data of the update variable set
  RealVector m_pData;

  /// unit vectors e_x and e_y
  std::vector<RealVector> m_unitVectors;

  /// fluxes at the plus and minus points of the stencil
  RealVector m_fluxPlus;
  RealVector m_fluxMinus;

  /// flux divergence with steps h and 2h
  RealVector m_divF;
  RealVector m_divF2;

  /// TNEQ source at the manufactured states of a cell
  RealVector m_srcAtMS;

  /// forcing per state and equation, indexed by local state ID
  std::vector<CFreal> m_forcing;

  /// true for the states whose forcing is computed
  std::vector<bool> m_isComputed;

  /// true once the checks of the first pass are logged
  bool m_checksDone;

  /// per equation: max |f|, max |S(U_MS)|, max |divF_h|, max |divF_h - divF_2h|
  RealVector m_maxForcing;
  RealVector m_maxSource;
  RealVector m_maxDivF;
  RealVector m_maxDivFDiff;

}; // class ManufacturedForcingTNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_ManufacturedForcingTNEQ_hh
