#ifndef COOLFluiD_SpectralFD_NullCommandSpectralFD_hh
#define COOLFluiD_SpectralFD_NullCommandSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that does nothing
 *
 * @author Kris Van den Abeele
 *
 */
class NullCommandSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit NullCommandSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~NullCommandSpectralFD();

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

}; // class NullCommandSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_NullCommandSpectralFD_hh
