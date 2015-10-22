#ifndef COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentReduce_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentReduce_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConcurrentCoupler/ConcurrentCouplerData.hh"
#include "ConcurrentCoupler/DataToTransfer.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a standard ConcurrentCoupler command that performs a 
 * global reduction of data
 *
 * @author Andrea Lani
 *
 */
class StdConcurrentReduce : public ConcurrentCouplerCom {
public:
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor.
   */
  explicit StdConcurrentReduce(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdConcurrentReduce();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Setup private data
   */
  virtual void setup();
  
  /**
   * Execute Processing actions
   */
  virtual void execute();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
protected:
  
  /// add information on the data to transfer
  /// @param idx  ID of the data socket to transfer
  void addDataToTransfer(const CFuint idx);
  
  /// create the group of processes for which reduction will be applied
  /// @param idx  ID of the data socket to transfer
  void createTransferGroup(const CFuint idx);
  
  /// reduce the data within the corresponding group of processes
  /// @param idx  ID of the data socket to transfer
  void reduceData(const CFuint idx);
  
protected: // data
  
  /// flag telling that the groups have been created
  bool m_createGroup;
  
  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;
  
  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// mapping from global to local IDs
  Common::CFMap<std::string, DataToTrasfer*> m_socketName2data;
  
  /// mapping from global to local IDs
  Common::CFMap<std::string, MPI_Op> m_nameToMPIOp;
    
  /// vector storing flags to identify ranks involved in the data transfer
  std::vector<std::vector<int> > m_isTransferRank;
  
  /// array of data to receive from other processes
  std::vector<CFreal> m_recvBuf;
  
  /// namespace involved in the reduction to which the current rank belongs
  std::vector<std::string> m_reduceNsp;
  
  /// names of the sending and receiving sockets with format:  
  /// "Namespace1_from>Namespace2_to" (no space on both sides of \">\".
  std::vector<std::string> m_socketsSendRecv;
  
  /// connectivity type for sockets (this is needed for defining global IDs)
  std::vector<std::string> m_socketsConnType;
  
  /// name of operation to be performed on each socket during reduction
  std::vector<std::string> m_operation;
    
}; // class StdConcurrentReduce
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentReduce_hh

