#ifndef COOLFluiD_Numerics_SubSystemCoupler_StdPreProcessWrite_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StdPreProcessWrite_hh

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
class StdPreProcessWrite : public CouplerCom {
public:

  /**
   * Constructor.
   */
  explicit StdPreProcessWrite(const std::string& name);

  /**
   * Destructor.
   */
  ~StdPreProcessWrite();

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

  void executeWriteOnTrs(const CFuint iProc);

  /**
   * In case the type of values to be passed is Nodal,
   * Create a datahandle with a vector of pointers to
   * the Geometric entities containing that node
   * @param socketName name of the socket that contains the NodeToFace connectivity
   */
  virtual void setNodeToFaceConnectivity(const std::string& socketName);

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the nodes
   * @param socketName name of the socket to store the coordinates of the nodes
   */
  virtual void fillNodalCoordDataHandle(const std::string& socketName);

  /**
   * In case the type of values to be passed is States,
   * Create a datahandle with a vector of pointers to
   * the Geometric entities containing that state
   * @param socketName name of the socket that contains the NodeToFace connectivity
   */
  virtual void setStateToFaceConnectivity(const std::string& socketName);

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the states
   * @param socketName name of the socket to store the coordinates of the nodes
   */
  virtual void fillStatesCoordDataHandle(const std::string& socketName);

private: // functions

  /**
   * Create the datahandle of coordinates and data
   * in case the transfer is done at the quadrature points
   */
  void setQuadratureDataHandles();

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the quadrature points
   * @param socketName name of the socket to store the coordinates of the quadrature points
   */
  void fillQuadratureCoordDataHandle(const std::string& socketName);

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the ghost points
   * @param socketName name of the socket to store the coordinates of the quadrature points
   */
  virtual void fillGhostCoordDataHandle(const std::string& socketName)
  {
    CFout << "StdPreProcessWrite::fillGhostCoordDataHandle(" << socketName <<") - You should not be here\n";
    cf_assert(false);
  }

protected: //

  void writeFile(const std::string& socketName);

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

  /// number of states coupled to another interface
  CFuint _nbCoupledStates;

  /// list of coordinate sockets already read
  std::vector<std::string> _alreadyReadSockets;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

}; // class StdPreProcessWrite

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StdPreProcessWrite_hh

