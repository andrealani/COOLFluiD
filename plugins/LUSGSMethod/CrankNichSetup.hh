#ifndef COOLFluiD_Numerics_LUSGSMethod_CrankNichSetup_hh
#define COOLFluiD_Numerics_LUSGSMethod_CrankNichSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "LUSGSMethod/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to setup the MeshData.
 *
 * @author Kris Van den Abeele
 */
class CrankNichSetup : public StdSetup {
public:

  /**
   * Constructor.
   */
  explicit CrankNichSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~CrankNichSetup();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// socket for rhs of past time step
  Framework::DataSocketSource< CFreal > socket_pastRhs;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_CrankNichSetup_hh
