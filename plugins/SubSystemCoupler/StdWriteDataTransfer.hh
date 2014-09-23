#ifndef COOLFluiD_Numerics_SubSystemCoupler_StdWriteDataTransfer_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StdWriteDataTransfer_hh

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
 * This class represents a StdWriteDataTransfer Coupler command
 *
 * @author Thomas Wuilbaut
 *
 */
class StdWriteDataTransfer : public CouplerCom {
public:

  /**
   * Constructor.
   */
  explicit StdWriteDataTransfer(const std::string& name);

  /**
   * Destructor.
   */
  ~StdWriteDataTransfer();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected:

  /**
   * Concrete execution of the command for each coupled processor.
   */
  virtual void executeWrite(const CFuint iProc);

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  ///Write a datahandle into a file
  void writeFile(const std::string filename);

protected: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  // handle to past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

}; // class StdWriteDataTransfer

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StdWriteDataTransfer_hh

