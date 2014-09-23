#ifndef COOLFluiD_SpectralFV_NullCommandSpectralFV_hh
#define COOLFluiD_SpectralFV_NullCommandSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that does nothing
 *
 * @author Kris Van den Abeele
 *
 */
class NullCommandSpectralFV : public SpectralFVMethodCom {
public:

  /**
   * Constructor.
   */
  explicit NullCommandSpectralFV(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~NullCommandSpectralFV();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup() {}

  /**
   * Unsetup private data
   */
  virtual void unsetup() {}

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args ) {}

  /**
   * Execute Processing actions
   */
  virtual void execute() {}

}; // class NullCommandSpectralFV

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_NullCommandSpectralFV_hh
