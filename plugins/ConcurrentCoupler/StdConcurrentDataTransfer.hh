#ifndef COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentDataTransfer_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentDataTransfer_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConcurrentCouplerData.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a standard ConcurrentCoupler command
 *
 * @author Andrea Lani
 *
 */
class StdConcurrentDataTransfer : public ConcurrentCouplerCom {
public:
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor.
   */
  explicit StdConcurrentDataTransfer(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdConcurrentDataTransfer();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // functions
  
  /// @return the DataStorage corresponding to the given namespace
  Common::SafePtr<Framework::DataStorage> getDataStorage(const std::string& nspName);
  
  /// gather data from all processes in namespace nspSend to 1 process in namespace nspRecv
  /// @param nspSend       namespace from which data are sent (>= 1 rank)
  /// @param nspRecv       namespace from which data are received (1 rank)
  /// @param sendSocketStr name of the socket from which data are sent (distributed)
  /// @param recvSocketStr name of the socket from which data are received (serial)
  void gatherData(const std::string& nspSend, 
		  const std::string& nspRecv,
		  const std::string& sendSocketStr, 
		  const std::string& recvSocketStr);
  
protected: // data
  
  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;
  
  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// names of the sending and receiving sockets with format:  
  /// "Namespace1_from>Namespace2_to" (no space on both sides of \">\".
  std::vector<std::string> _socketsSendRecv;
  
}; // class StdConcurrentDataTransfer
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentDataTransfer_hh

