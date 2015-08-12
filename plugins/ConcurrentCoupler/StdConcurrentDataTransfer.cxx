#include <numeric>

#include "Common/NotImplementedException.hh"
#include "Common/CFPrintContainer.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodCommand.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "ConcurrentCoupler/ConcurrentCouplerData.hh"
#include "ConcurrentCoupler/ConcurrentCoupler.hh"
#include "ConcurrentCoupler/StdConcurrentDataTransfer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdConcurrentDataTransfer, 
		      ConcurrentCouplerData, 
		      ConcurrentCouplerModule> 
stdConcurrentDataTransferProvider("StdConcurrentDataTransfer");

//////////////////////////////////////////////////////////////////////////////

void StdConcurrentDataTransfer::defineConfigOptions(Config::OptionList& options)
{  
  options.addConfigOption< vector<string> >
    ("SocketsSendRecv","Sockets to transfer, for example: Namespace1_send>Namespace2_recv (no space on both sides of \">\".");
}
  
//////////////////////////////////////////////////////////////////////////////

StdConcurrentDataTransfer::StdConcurrentDataTransfer(const std::string& name) :
  ConcurrentCouplerCom(name),
  _sockets(),
  socket_states("states")
{
  addConfigOptionsTo(this);
  
  _socketsSendRecv = vector<string>();
  setParameter("SocketsSendRecv",&_socketsSendRecv);
}
      
//////////////////////////////////////////////////////////////////////////////

StdConcurrentDataTransfer::~StdConcurrentDataTransfer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdConcurrentDataTransfer::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdConcurrentDataTransfer::configure ( Config::ConfigArgs& args )
{
  ConcurrentCouplerCom::configure(args);

  /*typedef std::vector<Common::SafePtr<CommandGroup> > InterfaceList;
  InterfaceList interfaces = getMethodData().getInterfaces();

  try {
    const std::string nsp = getMethodData().getNamespace();
  InterfaceList::iterator itr = interfaces.begin();
  for(; itr != interfaces.end(); ++itr) {

    const std::vector<std::string>& comNames = (*itr)->getComNames();

    // check if this command applies to this Interface
    if(count(comNames.begin(),comNames.end(),getName()) != 0) {

      const std::string interfaceName = (*itr)->getName();

      for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++iProc)
      {
        const vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);
        for (CFuint i=0; i< otherTrsNames.size(); ++i)
        {
          /// Get the datahandle for the data to be transfered TO THE OTHER SubSystem
          const std::string otherTrsName = otherTrsNames[i];

          vector<std::string> socketDataNames = getMethodData().getOtherCoupledDataName(interfaceName,otherTrsName,iProc);
          for(CFuint iType=0;iType < socketDataNames.size();iType++)
          {
            _sockets.createSocketSink<RealVector>(socketDataNames[iType]);
          } // coord types
        } // other trs
      } //loop over processors
    } // if check
  } // interfaces

  } catch (Exception& e)
  {
    CFout << e.what() << "\n" << CFendl;
    throw;
  }*/
}

//////////////////////////////////////////////////////////////////////////////

void StdConcurrentDataTransfer::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => start\n");
  
  const string couplingNsp = getMethodData().getNamespace();
  Group& group = PE::GetPE().getGroup(couplingNsp);
  // rank in MPI_COMM_WORLD
  const int rank  = PE::GetPE().GetRank("Default");
  // rank in group
  const int grank = PE::GetPE().GetRank(couplingNsp);
  const CFuint nbRanks = group.globalRanks.size();
  
  // number of variables that count for the coupling
  const CFuint couplingStride = PhysicalModelStack::getActive()->getNbEq();  
  
  SafePtr<VarSetTransformer> sendToRecvTrans = 
    getMethodData().getSendToRecvVecTrans();
  
  for (CFuint i = 0; i < _socketsSendRecv.size(); ++i) {
    vector<string> sendRecv = StringOps::getWords(_socketsSendRecv[i],'>');
    cf_assert(sendRecv.size() == 2);
    
    // namespace_socket (send)
    const string sendStr = sendRecv[0];
    vector<string> nspSocketSend = StringOps::getWords(sendStr,'_');
    cf_assert(nspSocketSend.size() == 2);
    
    // namespace_socket (recv)
    const string recvStr = sendRecv[1];
    vector<string> nspSocketRecv = StringOps::getWords(recvStr,'_');
    cf_assert(nspSocketRecv.size() == 2);
    
    const string nspSend    = nspSocketSend[0];
    const string socketSend = nspSocketSend[1];
    CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => send from " << nspSend << "-" << socketSend << "\n");
    const string nspRecv    = nspSocketRecv[0];
    const string socketRecv = nspSocketRecv[1];
    CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => recv from " << nspRecv << "-" << socketRecv << "\n");
    Group& groupSend = PE::GetPE().getGroup(nspSend);
    Group& groupRecv = PE::GetPE().getGroup(nspRecv);
    
    const CFuint nbRanksSend = groupSend.globalRanks.size();
    const CFuint nbRanksRecv = groupRecv.globalRanks.size();
        
    if (nbRanksSend > 1 && nbRanksRecv == 1) {
      CFuint sendcount = 0;
      vector<CFreal> recvbuf;
      vector<CFreal> sendbuf;
      vector<CFuint> sendIDs;
      vector<CFuint> recvIDs;
      vector<int> recvcounts(nbRanks, 0);
      vector<int> sendcounts(nbRanks, 0);
      vector<int> displs(nbRanks, 0);
      
      // this case gathers contributions from all ranks in the "send" namespace 
      // to a single rank correspoding to the "recv" namespace
      if (PE::GetPE().isRankInGroup(rank, nspSend)) {  
	recvbuf.resize(1);    // dummy in sending ranks
		
	// if my rank belong to the sending socket, gather first the number of elements to send
	SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
	  (SubSystemStatusStack::getCurrentName()).getNamespace(nspSend);
	SafePtr<MeshData> meshData = MeshDataStack::getInstance().getEntryByNamespace(nsp);
	SafePtr<DataStorage> ds = meshData->getDataStorage();
	CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => checking socket " << sendStr << "\n");
	const string sendLocal  = sendStr + "_local";
	const string sendGlobal = sendStr + "_global";
	
	// local data (CFreal)
	if (ds->checkData(sendStr)) {
	  DataHandle<CFreal> array = ds->getData<CFreal>(sendStr);
	  CFLog(VERBOSE, "P" << rank << " has socket " << sendStr << "\n"); 
	}
	// global data (State*)
	else if (ds->checkData(sendLocal) && ds->checkData(sendGlobal)) {
	  DataHandle<State*, GLOBAL> array = ds->getGlobalData<State*>(sendStr);
	  CFLog(VERBOSE, "P" << rank << " has socket " << sendStr << " with sizes = [" 
		<< array.getLocalSize() << ", " << array.getGlobalSize() << "]\n"); 
	  sendcount = array.getLocalSize()*couplingStride;
	  // sendbuf = array.getGlobalData(0);
	  // or = &array.getGlobalData(0)[0]; ??
	  sendbuf.reserve(sendcount);
	  sendIDs.reserve(sendcount);
	  
	  for (CFuint ia = 0; ia < array.size(); ++ia) {
	    const CFuint globalID = array[ia]->getGlobalID();
	    if (array[ia]->isParUpdatable()) {
	      const State *const tState = sendToRecvTrans->transform(array[ia]);
	      cf_assert(tState->size() == couplingStride);
	      
	      for (CFuint s = 0; s < couplingStride; ++s) {
		sendbuf.push_back((*tState)[s]);
		sendIDs.push_back(globalID*couplingStride+s);
	      }
	    }
	  }
	  cf_assert(sendbuf.size() == sendcount);
	}
	
	// fill in the number of counts to send from this rank
	sendcounts[grank] = sendcount;
      }
            
      MPIError::getInstance().check
	("MPI_Allreduce", "StdConcurrentDataTransfer::execute()", 
	 MPI_Allreduce(&sendcounts[0], &recvcounts[0], nbRanks,
		       MPIStructDef::getMPIType(&recvcounts[0]), MPI_MAX, group.comm));
      
      CFLog(DEBUG_MAX, CFPrintContainer<vector<int> >
	    ("StdConcurrentDataTransfer::execute() => recvcounts  = ", &recvcounts));
      
      displs[0] = 0;
      CFuint totRecvcount = recvcounts[0];
      for (CFuint r = 1; r < nbRanks; ++r) {
	if (recvcounts[r] > 0) {
	  displs[r] = totRecvcount;
	}
	totRecvcount += recvcounts[r];
      }
      cf_assert(totRecvcount == std::accumulate(recvcounts.begin(), recvcounts.end(),0));
      
      int root = -1;
      int sendroot = -1; 
      if (PE::GetPE().isRankInGroup(rank, nspRecv)) {
	sendroot = PE::GetPE().GetRank(couplingNsp);
	recvbuf.resize(totRecvcount);
	sendIDs.resize(1);
	recvIDs.resize(totRecvcount);
      }
      
      MPIError::getInstance().check
	("MPI_Allreduce", "StdConcurrentDataTransfer::execute()", 
	 MPI_Allreduce(&sendroot, &root, 1, MPIStructDef::getMPIType(&root), MPI_MAX, group.comm));
      cf_assert(root >= 0);
      
      // MPIStruct sendMS;
      // int lns[2];
      // lns[0] = lns[1] = sendcount;
      // MPIStructDef::buildMPIStruct<CFuint,CFreal>(&sendIDs[0], &sendbuf[0], lns, sendMS);
      
      // MPIStruct recvMS;
      // int lnr[2];
      // lnr[0] = lnr[1] = totRecvcount;
      // MPIStructDef::buildMPIStruct<CFuint,CFreal>(&recvIDs[0], &recvbuf[0], lnr, recvMS);
      
      // MPIError::getInstance().check
      // 	("MPI_Gatherv", "StdConcurrentDataTransfer::execute()", 
      // 	 MPI_Gatherv(sendMS.start, 1, sendMS.type, recvMS.start, &recvcounts[0], &displs[0], 
      // 		     recvMS.type, root, group.comm));
      
      MPIError::getInstance().check
      	("MPI_Gatherv", "StdConcurrentDataTransfer::execute()", 
      	 MPI_Gatherv(&sendbuf[0], sendcount, MPIStructDef::getMPIType(&sendbuf[0]),
      		     &recvbuf[0], &recvcounts[0], &displs[0], 
      		     MPIStructDef::getMPIType(&sendbuf[0]), root, group.comm));
      
      MPIError::getInstance().check
      	("MPI_Gatherv", "StdConcurrentDataTransfer::execute()", 
      	 MPI_Gatherv(&sendIDs[0], sendcount, MPIStructDef::getMPIType(&sendIDs[0]),
      		     &recvIDs[0], &recvcounts[0], &displs[0], 
      		     MPIStructDef::getMPIType(&sendIDs[0]), root, group.comm));
      
      if (grank == root) {
	// order the received data
	SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
	  (SubSystemStatusStack::getCurrentName()).getNamespace(nspRecv);
	SafePtr<MeshData> meshData = MeshDataStack::getInstance().getEntryByNamespace(nsp);
	SafePtr<DataStorage> ds = meshData->getDataStorage();
	CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => checking socket " << recvStr << "\n");
	const string recvLocal  = recvStr + "_local";
	const string recvGlobal = recvStr + "_global";
	
	// local data (CFreal)
	if (ds->checkData(recvStr)) {
	  DataHandle<CFreal> array = ds->getData<CFreal>(recvStr);
	  CFLog(VERBOSE, "P" << rank << " has socket " << recvStr << "\n"); 
	}
	// global data (State*)
	else if (ds->checkData(recvLocal) && ds->checkData(recvGlobal)) {
	  DataHandle<State*, GLOBAL> array = ds->getGlobalData<State*>(recvStr);
	  CFLog(VERBOSE, "P" << rank << " has socket " << recvStr << " with sizes = [" 
		<< array.getLocalSize() << ", " << array.getGlobalSize() << "]\n"); 
	  cf_assert(array.getLocalSize() == array.getGlobalSize());
	  cf_assert(array.size() == array.getGlobalSize());
	  cf_assert(array.size()*couplingStride == totRecvcount);
	  
	  CFreal* sarray = array.getGlobalArray()->ptr();
	  for (CFuint is = 0; is < totRecvcount; ++is) {
	    cf_assert(is < recvIDs.size());
	    cf_assert(is < recvbuf.size());
	    sarray[recvIDs[is]] = recvbuf[is];  
	  }
	}
      }
    }
    else if (nbRanksSend == 1 && nbRanksRecv > 1) {
      // use MPI_Scatterv
    }
    else if (nbRanksSend > 1 && nbRanksRecv > 1) {
      // use MPI_Alltoallv
    }
  }
  
  CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

  } // namespace Numerics

} // namespace COOLFluiD
