#ifndef COOLFluiD_Numerics_SubSystemCoupler_StdPreProcessRead_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StdPreProcessRead_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Framework/DynamicDataSocketSet.hh"
#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to setup the SubSystemCoupler Method
 *
 * @author Thomas Wuilbaut
 *
 */
class StdPreProcessRead : public CouplerCom {
public:

  /**
   * Constructor.
   */
  explicit StdPreProcessRead(const std::string& name);

  /**
   * Destructor.
   */
  ~StdPreProcessRead();

  /**
   * Configures this command
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute Processing actions
   */
  void executeOnTrs();

  /**
   * Execute Processing actions
   */
  void executeReadOnTrs(const CFuint iProc);

protected: //

  //Read file containing the coordinates of the other side
  void readFile(const std::string& name);

  /**
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);

protected: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// @todo missing documentation
  CFuint _nbOtherStates;

  /// list of coordinate sockets already read
  std::vector<std::string> _alreadyReadSockets;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for State's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

}; // class StdPreProcessRead

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StdPreProcessRead_hh

