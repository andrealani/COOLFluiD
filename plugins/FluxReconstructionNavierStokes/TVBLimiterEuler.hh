#ifndef COOLFluiD_FluxReconstructionMethod_TVBLimiterEuler_hh
#define COOLFluiD_FluxReconstructionMethod_TVBLimiterEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/TVBLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that applies an elementwise TVB limiter to the solution,
 * taking into account the requirement of pressure positivty (pressure is limited instead of rhoE)
 *
 * @author Kris Van den Abeele
 *
 */
class TVBLimiterEuler : public TVBLimiter {
public:

  /**
   * Constructor.
   */
  explicit TVBLimiterEuler(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~TVBLimiterEuler();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

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
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // functions

  /**
   * compute the cell averaged state
   */
  virtual void reconstructCellAveragedState();

  /**
   * compute a cell averaged variable
   */
  virtual void reconstructCellAveragedVariable(const CFuint iEq);

  /**
   * compute derivative in the cell center of a variable
   */
  virtual void computeCellCenterDerivVariable(const CFuint iEq);

  /**
   * check if limiting is necessary
   * @pre solution has been reconstructed at flux points
   */
  virtual void setLimitBooleans();

  /**
   * apply the limiter
   */
  virtual void limitStates();

protected: // functions

  /// compute the pressure from the conserved variables
  virtual void computePressFromConsVar(const RealVector& consVar, CFreal& press) = 0;

  /// compute the rhoE from pressure and the other conserved variables
  virtual void computeRhoEFromPressAndOtherConsVar(RealVector& consVar, const CFreal& press) = 0;

protected: // data

  /// number of equations minus one
  CFuint m_nbrEqsMinusOne;

  /// boolean telling whether to check for density/pressure positivity
  bool m_positivityCheck;

  /// minimum allowable value for density
  CFreal m_minDensity;

  /// minimum allowable value for pressure
  CFreal m_minPressure;

}; // class TVBLimiterEuler

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_TVBLimiterEuler_hh
