#ifndef COOLFluiD_SpectralFV_DivideByVolFracRHSSpectralFV_hh
#define COOLFluiD_SpectralFV_DivideByVolFracRHSSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFV/BaseVolTermComputer.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that multiplies the residuals with the inverse CV volume fractions
 *
 * @author Kris Van den Abeele
 *
 */
class DivideByVolFracRHSSpectralFV : public SpectralFVMethodCom {
public:

  /**
   * Constructor.
   */
  explicit DivideByVolFracRHSSpectralFV(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~DivideByVolFracRHSSpectralFV();

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

  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /**
   * Divides the rhs by the cell volume
   */
  void divideRHSByVolume();

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// inverse volume fractions of CVs
  Common::SafePtr< std::vector< CFreal > > m_invVolFracCVs;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class DivideByVolFracRHSSpectralFV

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_DivideByVolFracRHSSpectralFV_hh
