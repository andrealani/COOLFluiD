#ifndef COOLFluiD_Numerics_LUSGSMethod_BDF3Setup_hh
#define COOLFluiD_Numerics_LUSGSMethod_BDF3Setup_hh

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
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class BDF3Setup : public StdSetup {
public:

  /**
   * Constructor.
   */
  explicit BDF3Setup(const std::string& name);

  /**
   * Destructor.
   */
  ~BDF3Setup();

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

  /// socket for the past time rhs
  Framework::DataSocketSource< CFreal > socket_pastTimeRhs;

  /// socket for the past past time rhs
  Framework::DataSocketSource< CFreal > socket_pastPastTimeRhs;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LUSGSMethod_BDF3Setup_hh
