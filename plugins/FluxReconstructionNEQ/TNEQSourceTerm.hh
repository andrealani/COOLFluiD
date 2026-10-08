#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_TNEQSourceTerm_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_TNEQSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNEQ/CNEQSourceTerm.hh"

#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * A command for adding the source term for TNEQ
 *
 * @author Ray Vandenhoeck
 * @author Firas Ben Ameur
 * @author Rayan Dhib
 *
 */
class TNEQSourceTerm : public CNEQSourceTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit TNEQSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~TNEQSourceTerm();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * get data required for source term computation
   */
  void getSourceTermData();

  /**
   * add the source term
   */
  void addSourceTerm(RealVector& resUpdates);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set the adimensional vibrational temperatures
   */
  virtual void setVibTemperature(const RealVector& pdata, 
				 const Framework::State& state,
				 RealVector& tvib);
      
  /**
   * Electron pressure work at one solution point (option ElectronPressureWork): sets m_pe = rho_e R_e T_e
   * [Pa] and m_divV = div(u) [1/s] (plus v/r with Data.Axisymmetric), with rho_e = rho y_e, R_e = R_u/M_e
   * and T_e the temperature of the equation that carries the free-electron energy
   * @param iSol index of the solution point in the current cell
   * @param rhodim dimensional density at that point
   */
  void computePeDivV(const CFuint iSol, const CFreal rhodim);

  /// Compute the energy transfer terms
  void computeSourceVT(RealVector& omegaTv, CFreal& omegaRad);

  /**
   * Analytical Jacobian of the source term at one solution point (AnalyticalJacob = true).
   * The library gives J(i,j) = d prodterm_i / d W_j, with W = [rho_s, momentum, T, Tv] in SI units.
   * The chain rule to the update variables U is diagonal: dW_j/dU_j = rho_ref, T_ref (RhoivtTv)
   * or rho_s, T, Tv (LogRhoivLogTTv, U = ln W).
   * Stores m_stateJacobian[col][row] = dR_row/dU_col, R being the residual update of addSourceTerm.
   * @param iState index of the solution point in the current cell
   */
  virtual void getSToStateJacobian(const CFuint iState);

  /// true if the source term has to be computed for the active equation subsystem
  bool doComputeSourceTerm() const;

  /// fill the dimensional inputs of the chemistry library (T, Tv, p, rho, y) at one solution point
  void setLibraryInputs(const Framework::State& state, CFreal& pdim, CFreal& Tdim, CFreal& rhodim);

protected: // data

  /// radiation source term
  CFdouble m_omegaRad;
  
  /// velocity divergence
  CFreal m_divV;
  
  /// electron pressure
  CFreal m_pe;
  
  /// array to store the energy exchange terms
  RealVector m_omegaTv;
  
  /// reference data
  Common::SafePtr<RealVector> m_refData;

  /// true if the update state holds logarithmic variables (LogRhoivLogTTv: ln rho_i, u, v, ln T, ln Tv)
  bool m_logVariables;

  /// PLATO source Jacobian d prodterm_i / d W_j, W = [rho_s, momentum, T, Tv] (SI)
  RealMatrix m_platoJacob;

  /// option: add the electron pressure work -p_e div(u) to the free-electron energy equation
  bool m_electronPressureWork;

  /// true if the electron pressure work is computed (option on and free electrons in the mixture)
  bool m_addPeDivV;

  /// index in the temperature vector m_tvDim of the free-electron temperature T_e
  CFuint m_teID;

  /// gas constant of the free electrons R_e = R_u / M_e [J/(kg K)]
  CFreal m_Re;

  /// socket for the gradients of the update variables at the solution points
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

}; // class TNEQSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_TNEQSourceTerm_hh

