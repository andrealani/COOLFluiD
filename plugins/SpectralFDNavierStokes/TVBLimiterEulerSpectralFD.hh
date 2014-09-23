#ifndef COOLFluiD_SpectralFD_TVBLimiterEulerSpectralFD_hh
#define COOLFluiD_SpectralFD_TVBLimiterEulerSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/ReconstructStatesSpectralFD.hh"
#include "SpectralFD/SpectralFDMethodData.hh"
#include "SpectralFD/TVBLimiterSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that applies an elementwise TVB limiter to the solution,
 * taking into account the requirement of pressure positivty (pressure is limited instead of rhoE)
 *
 * @author Kris Van den Abeele
 *
 */
class TVBLimiterEulerSpectralFD : public TVBLimiterSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit TVBLimiterEulerSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~TVBLimiterEulerSpectralFD();

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

}; // class TVBLimiterEulerSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_TVBLimiterEulerSpectralFD_hh
