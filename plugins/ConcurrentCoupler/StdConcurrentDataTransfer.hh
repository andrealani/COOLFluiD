#ifndef COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentDataTransfer_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentDataTransfer_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConcurrentCouplerData.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class VarSetTransformer;
  }
  
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
  
protected: // functions
  
  /// struct holding some data for controlling the transfer
  class DataToTrasfer {
  public:
    DataToTrasfer() 
    {array = CFNULL; sendStride = recvStride = nbRanksSend = nbRanksRecv = 0;}
    
    ~DataToTrasfer() {}
    
    CFreal* array;        // local array (send or recv)
    CFuint arraySize;     // total size of the local array<CFreal> (send or recv)
    CFuint sendStride;    // stride for the send socket
    CFuint recvStride;    // stride for the recv socket  
    CFuint nbRanksSend;   // number of ranks in the send group   
    CFuint nbRanksRecv;   // number of ranks in the recv group
    std::string dofsName; // name of the corresponding dofs ("*_states" or "*_nodes")
    std::string nspSend;      // namespace from which data are sent
    std::string nspRecv;      // namespace from which data are received
    std::string sendSocketStr; // name of the socket from which data are sent
    std::string recvSocketStr; // name of the socket from which data are received 
  };
  
  /// @return the DataStorage corresponding to the given namespace
  Common::SafePtr<Framework::DataStorage> getDataStorage(const std::string& nspName);
  
  /// gather data from all processes in namespace nspSend to 1 process in namespace nspRecv
  /// @param idx           index of the data transfer
  virtual void gatherData(const CFuint idx);
  
  /// scatter data from 1 process in namespace nspSend to all processes in namespace nspRecv
  /// @param idx           index of the data transfer
  virtual void scatterData(const CFuint idx);
  
  /// fill a mapping between global and local IDs
  /// @param ds            pointer to DataStorage
  /// @param socketName    name of the socket
  /// @param global2local  reference to the map object to fill in
  template <typename T>
  void fillMapGlobalToLocal(Common::SafePtr<DataToTrasfer> dtt,
			    Common::SafePtr<Framework::DataStorage> ds,
			    Common::CFMap<CFuint, CFuint>& global2local);
  
  /// fill data to be sent
  /// @param vsTrans         variable transformer
  /// @param ds              pointer to DataStorage
  /// @param dofStr          name of the dof socket
  /// @param sendcount       size of buffer to send
  /// @param sendbuf         buffer of data to be sent
  /// @param sendIDs         buffer od global IDs to be sent
  /// @param array           local array from which data will be sent
  /// @param arraySize       size of local array from which data will be sent
  /// @param isDof           tells if this is the dof array
  template <typename T>
  void fillSendDataGather(Common::SafePtr<DataToTrasfer> dtt,
			  Common::SafePtr<Framework::VarSetTransformer> vsTrans,
			  Common::SafePtr<Framework::DataStorage> ds,
			  CFuint& sendcount, 
			  std::vector<CFreal>& sendbuf,
			  std::vector<CFuint>& sendIDs);
  
  /// fill data to be sent
  /// @param ds              pointer to DataStorage
  /// @param dofStr          name of the dof socket
  /// @param sendcounts      size of buffer to send in each rank
  /// @param sendIDcounts    size of IDs to send in each rank
  /// @param array           local array from which data will be sent
  /// @param arraySize       size of local array from which data will be sent
  /// @param isDof           tells if this is the dof array
  template <typename T>
  void fillSendCountsScatter(Common::SafePtr<DataToTrasfer> dtt,
			     Common::SafePtr<Framework::DataStorage> ds,
			     std::vector<int>& sendcounts,
			     std::vector<int>& sendIDcounts); 
  
  /// @return the rank (within nspCoupling) of the root process belonging to namespace nsp
  /// @param nsp           namespace to which the process belongs
  /// @param nspCoupling   coupling namespace 
  int getRootProcess(const std::string& nsp, const std::string& nspCoupling) const;
  
  /// add information on the data to transfer
  /// @param idx  ID of the data socket to transfer
  void addDataToTransfer(const CFuint idx);
  
protected: // data
  
  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;
  
  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// vector transformer from send (source) to recv (target) variables
  std::vector<Common::SelfRegistPtr<Framework::VarSetTransformer> > _sendToRecvVecTrans;
  
  /// mapping from global to local IDs
  Common::CFMap<CFuint, CFuint> _global2localIDs;
  
  /// mapping from global to local IDs
  Common::CFMap<std::string, DataToTrasfer*> _socketName2data;
  
  /// names of the sending and receiving sockets with format:  
  /// "Namespace1_from>Namespace2_to" (no space on both sides of \">\".
  std::vector<std::string> _socketsSendRecv;
  
  /// connectivity type for sockets (this is needed for defining global IDs)
  std::vector<std::string> _socketsConnType;
  
  /// variables transformers from send to recv variables
  std::vector<std::string> _sendToRecvVecTransStr;
  
}; // class StdConcurrentDataTransfer
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ConcurrentCoupler/StdConcurrentDataTransfer.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ConcurrentCoupler_StdConcurrentDataTransfer_hh

