#ifndef COOLFluiD_Numerics_ConcurrentCoupler_DataToTransfer_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_DataToTransfer_hh

//////////////////////////////////////////////////////////////////////////////

#include <mpi.h>

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a tuple of data to define the transfer during 
 * parallel communication for a concurrent coupler
 *
 * @author Andrea Lani
 */
class DataToTrasfer {
public:
  
  /// default constructor
  DataToTrasfer() 
  {array = CFNULL; sendStride = recvStride = nbRanksSend = nbRanksRecv = 0;}
  
  // destructor
  ~DataToTrasfer() {}
  
  CFreal* array;             // local array (send or recv)
  CFuint arraySize;          // total size of the local array<CFreal> (send or recv)
  CFuint sendStride;         // stride for the send socket
  CFuint recvStride;         // stride for the recv socket  
  CFuint nbRanksSend;        // number of ranks in the send group   
  CFuint nbRanksRecv;        // number of ranks in the recv group
  std::string dofsName;      // name of the corresponding dofs ("*_states" or "*_nodes")
  std::string nspSend;       // namespace from which data are sent
  std::string nspRecv;       // namespace from which data are received
  std::string sendSocketStr; // name of the socket from which data are sent
  std::string recvSocketStr; // name of the socket from which data are received 
  std::string groupName;     // name of the MPI group in which data transfer is active
  MPI_Op operation;          // MPI operation to apply
};
 
//////////////////////////////////////////////////////////////////////////////
    
    } // namespace ConcurrentCoupler
    
  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ConcurrentCoupler_DataToTransfer_hh
