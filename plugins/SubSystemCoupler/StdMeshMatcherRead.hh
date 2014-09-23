#ifndef COOLFluiD_Numerics_SubSystemCoupler_StdMeshMatcherRead_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StdMeshMatcherRead_hh

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
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to set the
   * match between meshes
   *
   * @author Thomas Wuilbaut
   *
   */

class StdMeshMatcherRead : public CouplerCom {
public:

  /**
   * Constructor.
   */
  explicit StdMeshMatcherRead(const std::string& name);

  /**
   * Destructor.
   */
  ~StdMeshMatcherRead();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Executes the command.
   */
  void execute();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data

  void executeRead();

  /**
   * Reading from a file the acceptance status of the points
   */
  virtual void readIsAcceptedFile(const std::string dataHandleName);

  /**
   * Assembling the acceptance status of the states from the various
   * processors in a unique datahandle
   */
  virtual void assembleIsAcceptedFiles(const std::string localDataHandleName);

  /**
   * Resize the coupled data datahandles
   */
  void resizeDataHandles();

  /**
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  virtual void getWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);

protected: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// Index of the processor being read
  CFuint _iProc;

  /// Temporary variable to count the number of accepted states
  RealVector _tempIsAccepted;

  /// Temporary variable to store which processor accepted the states
  std::vector<CFuint> _tempParallelIndex;


}; // class StdMeshMatcherRead

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StdMeshMatcherRead_hh

