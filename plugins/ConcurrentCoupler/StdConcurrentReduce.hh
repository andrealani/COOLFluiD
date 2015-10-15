#ifndef COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentReduce_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentReduce_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>

#include "ConcurrentCoupler/ConcurrentCouplerData.hh"
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
  void createReduceGroup();
  
  /// @return true if the current rank has to be involved in reduction
  bool isReductionRank()
  {
    const int nspRank  = 
      Common::PE::GetPE().GetRank(getMethodData().getNamespace());
    cf_assert(nspRank < m_isReduceRank.size());
    return (m_isReduceRank[nspRank] == 1); 
  }
  
protected:
  
  /// struct holding some data for controlling the transfer
  class DataToTrasfer {
  public:
    /// constructor
    DataToTrasfer() {array = CFNULL; sendStride = nbRanksSend = 0;}
    
    ///destructor
    ~DataToTrasfer() {}
    
    CFreal* array;             // local array (send or recv)
    CFuint arraySize;          // total size of the local array<CFreal> (send or recv)
    CFuint sendStride;         // stride for the send socket
    CFuint nbRanksSend;        // number of ranks in the send group   
    std::string dofsName;      // name of the corresponding dofs ("*_states" or "*_nodes")
    std::string nspSend;       // namespace from which data are sent
    std::string sendSocketStr; // name of the socket from which data are sent
    MPI_Op operation;          // MPI operation to apply
  };
  
protected: // data
  
  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;
  
  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// mapping from global to local IDs
  Common::CFMap<std::string, DataToTrasfer*> m_socketName2data;
  
  /// mapping from global to local IDs
  Common::CFMap<std::string, MPI_Op> m_nameToMPIOp;
  
  /// integer flag telling if this is a rank involved in the reduction
  std::vector<int> m_isReduceRank;
  
  /// namespace involved in the reduction to which te current rank belongs
  std::string m_reduceNsp;
  
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

