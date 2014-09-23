#ifndef COOLFluiD_Numerics_SubSystemCoupler_StdReadDataTransfer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StdReadDataTransfer_hh

//////////////////////////////////////////////////////////////////////////////

#include "SubSysCouplerData.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a StdReadDataTransfer Coupler command
 *
 * @author Thomas Wuilbaut
 *
 */
class StdReadDataTransfer : public CouplerCom {
public:

  /**
   * Constructor.
   */
  explicit StdReadDataTransfer(const std::string& name);

  /**
   * Destructor.
   */
  ~StdReadDataTransfer();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Sets up the command.
   */
  virtual void setup();

protected: // functions

  /**
   * Execute Processing actions
   */
  void executeRead();

  /**
   * Apply the post transformer on the data read
   */
  void transformReceivedData();

  /**
   * Read the data from a file and put the value in the data datahandle
   */
  void readFile(const std::string filename, const std::string acceptedFileName,const std::string dataHandleName,const std::string acceptedDataHandleName);

  ///Read the datahandle of the other namespace put values into your own datahandle
  void readFromDataHandle(const std::string dataHandleName);

  ///Outputs to file the norm of the data update
  void prepareNormFile(const std::string dataHandleName);

  ///Outputs to file the norm of the data update
  void computeNorm(const std::string dataHandleName);

  /**
   * Use the PostVariableTransformer to transform
   * the values received through the datahandle
   */
  virtual void transformReceivedGhostData()
  {
    CFout << "StdReadDataTransfer::transformReceivedGhostData - You should not be here\n";
    CFout << "use FVMCCReadDataTransfer\n";
    cf_assert(false);
  }

  /**
   * Use the PostVariableTransformer to transform
   * the values received through the datahandle
   */
  void transformReceivedQuadratureData();

  /**
   * Use the PostVariableTransformer to transform
   * the values received through the datahandle
   */
  virtual void transformReceivedNodalData();

  /**
   * Use the PostVariableTransformer to transform
   * the values received through the datahandle
   */
  void transformReceivedStatesData();

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

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  ///temp storage of interface name
  std::string _interfaceName;

  ///temp storage of currentTRS name
  std::string _currentTrsName;

  ///other processor for which the data is processed
  CFuint _iProc;

}; // class StdReadDataTransfer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StdReadDataTransfer_hh

