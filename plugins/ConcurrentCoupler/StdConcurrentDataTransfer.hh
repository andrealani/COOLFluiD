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

