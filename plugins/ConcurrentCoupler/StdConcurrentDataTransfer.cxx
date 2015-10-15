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
#include "Framework/VarSetTransformer.hh"

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
  options.addConfigOption< vector<string> >
    ("SocketsConnType","Connectivity type for sockets to transfer (State or Node): this is ne1eded to define global IDs.");
  options.addConfigOption< vector<string> >
    ("SendToRecvVariableTransformer","Variables transformers from send to recv variables.");
}
      
//////////////////////////////////////////////////////////////////////////////

StdConcurrentDataTransfer::StdConcurrentDataTransfer(const std::string& name) :
  ConcurrentCouplerCom(name),
  _sockets(),
  socket_states("states"),
  _sendToRecvVecTrans(),
  _global2localIDs(),
  _socketName2data()
{
  addConfigOptionsTo(this);
  
  _socketsSendRecv = vector<string>();
  setParameter("SocketsSendRecv",&_socketsSendRecv);
  
  _socketsConnType = vector<string>();
  setParameter("SocketsConnType",&_socketsConnType);
  
  _sendToRecvVecTransStr = vector<string>();
  setParameter("SendToRecvVariableTransformer", &_sendToRecvVecTransStr);
}
      
//////////////////////////////////////////////////////////////////////////////

StdConcurrentDataTransfer::~StdConcurrentDataTransfer()
{
  for (CFuint i =0 ; i < _socketName2data.size(); ++i) {
    delete _socketName2data[i];
  }  
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdConcurrentDataTransfer::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdConcurrentDataTransfer::configure ( Config::ConfigArgs& args )
{
  ConcurrentCouplerCom::configure(args);
  
  if (_socketsConnType.size() != _socketsSendRecv.size()) {
    CFLog(ERROR, "StdConcurrentDataTransfer::configure() => SocketsSendRecv.size() != SocketsConnType.size()\n");
    cf_assert(_socketsConnType.size() == _socketsSendRecv.size());
  }
  
  // configure variable transformers
  const string name = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = 
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(name);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  SafePtr<VarSetTransformer::PROVIDER> vecTransProv = CFNULL;
  
  if (_sendToRecvVecTransStr.size() == 0) {
    _sendToRecvVecTransStr.resize(_socketsSendRecv.size(), "Identity");
  }
  _sendToRecvVecTrans.resize(_sendToRecvVecTransStr.size());
  
  for (CFuint i = 0; i < _sendToRecvVecTransStr.size(); ++i) {
    CFLog(VERBOSE, "Configuring VarSet Transformer: " << _sendToRecvVecTransStr[i] << "\n");
    
    try {
      vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider
	(_sendToRecvVecTransStr[i]);
    }
    catch (Common::NoSuchValueException& e) {
      _sendToRecvVecTransStr[i] = "Identity";
      
      CFLog(VERBOSE, e.what() << "\n");
      CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
      vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider
	(_sendToRecvVecTransStr[i]);
    }
    
    cf_assert(vecTransProv.isNotNull());  
    _sendToRecvVecTrans[i].reset(vecTransProv->create(physModel->getImplementor()));
    cf_assert(_sendToRecvVecTrans[i].getPtr() != CFNULL);
  }
}

//////////////////////////////////////////////////////////////////////////////
  
void StdConcurrentDataTransfer::setup()
{
  // set up the variable transformer
  for (CFuint i = 0; i < _sendToRecvVecTrans.size(); ++i){
    _sendToRecvVecTrans[i]->setup(1);
  }
}
 
//////////////////////////////////////////////////////////////////////////////
       
void StdConcurrentDataTransfer::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => start\n");
  
  // this should go in the setup, but it uses blocking MPI collective calls 
  // here is less harmful 
  if (_socketName2data.size() == 0) {
    // create a preliminary mapping between sockets names and related data to transfer
    for (CFuint i = 0; i < _socketsSendRecv.size(); ++i) {
      addDataToTransfer(i);
    }
  }
  
  for (CFuint i = 0; i < _socketsSendRecv.size(); ++i) {
    SafePtr<DataToTrasfer> dtt = _socketName2data.find(_socketsSendRecv[i]); 
    const CFuint nbRanksSend = dtt->nbRanksSend;
    const CFuint nbRanksRecv = dtt->nbRanksRecv;
    
    if (nbRanksSend > 1 && nbRanksRecv == 1) {
      gatherData(i);
    }
    else if (nbRanksSend == 1 && nbRanksRecv > 1) {
      scatterData(i);
    }
    else if (nbRanksSend > 1 && nbRanksRecv > 1) {
      throw NotImplementedException
	(FromHere(),"StdConcurrentDataTransfer::execute() => (nbRanksSend > 1 && nbRanksRecv > 1)");
    }
  }
  
  CFLog(VERBOSE, "StdConcurrentDataTransfer::execute() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void StdConcurrentDataTransfer::gatherData(const CFuint idx)
{
  SafePtr<DataToTrasfer> dtt = _socketName2data.find(_socketsSendRecv[idx]); 
  
  const string nspSend = dtt->nspSend;
  const string nspRecv = dtt->nspRecv;
  const string sendSocketStr = dtt->sendSocketStr;
  const string recvSocketStr = dtt->recvSocketStr;
  const string nspCoupling = getMethodData().getNamespace();
  
  CFLog(INFO, "StdConcurrentDataTransfer::gatherData() from namespace[" << nspSend 
	<< "] to namespace [" << nspRecv << "] => start\n");
  
  Group& group    = PE::GetPE().getGroup(nspCoupling);
  const int rank  = PE::GetPE().GetRank("Default"); // rank in MPI_COMM_WORLD
  const int grank = PE::GetPE().GetRank(nspCoupling); // rank in group
  const CFuint nbRanks = group.globalRanks.size();
  
  // number of variables that count for the coupling
  cf_assert(idx < _sendToRecvVecTrans.size());
  SafePtr<VarSetTransformer> sendToRecvTrans = _sendToRecvVecTrans[idx].getPtr();
  cf_assert(sendToRecvTrans.isNotNull());
  
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
    recvbuf.resize(1); // dummy in sending ranks
		
    // if my rank belong to the sending socket, gather first the number of elements to send
    SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspSend);
    if (_socketsConnType[idx] == "State") {
      fillSendDataGather<State*>(dtt, sendToRecvTrans, ds, sendcount, sendbuf, sendIDs);
    }
    if (_socketsConnType[idx] == "Node") {
      fillSendDataGather<Node*>(dtt, sendToRecvTrans, ds, sendcount, sendbuf, sendIDs);
    }
    
    cf_assert(sendbuf.size() == sendcount);
    
    // fill in the number of counts to send from this rank
    sendcounts[grank] = sendcount;
  }
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "StdConcurrentDataTransfer::gatherData()", 
     MPI_Allreduce(&sendcounts[0], &recvcounts[0], nbRanks,
		   MPIStructDef::getMPIType(&recvcounts[0]), MPI_MAX, group.comm));
  
  CFLog(DEBUG_MAX, CFPrintContainer<vector<int> >
	("StdConcurrentDataTransfer::gatherData() => recvcounts  = ", &recvcounts));
  
  displs[0] = 0;
  CFuint totRecvcount = recvcounts[0];
  for (CFuint r = 1; r < nbRanks; ++r) {
    if (recvcounts[r] > 0) {
      displs[r] = totRecvcount;
    }
    totRecvcount += recvcounts[r];
  }
  cf_assert(totRecvcount == std::accumulate(recvcounts.begin(), recvcounts.end(),0));
  
  if (PE::GetPE().isRankInGroup(rank, nspRecv)) {
    recvbuf.resize(totRecvcount);
    sendIDs.resize(1);
    recvIDs.resize(totRecvcount);
  }
  
  int root = getRootProcess(nspRecv, nspCoupling);
  
  // transfer the actual data
  MPIError::getInstance().check
    ("MPI_Gatherv", "StdConcurrentDataTransfer::gatherData()", 
     MPI_Gatherv(&sendbuf[0], sendcount, MPIStructDef::getMPIType(&sendbuf[0]),
		 &recvbuf[0], &recvcounts[0], &displs[0], 
		 MPIStructDef::getMPIType(&sendbuf[0]), root, group.comm));
  
  // transfer the global IDs
  MPIError::getInstance().check
  ("MPI_Gatherv", "StdConcurrentDataTransfer::gatherData()", 
   MPI_Gatherv(&sendIDs[0], sendcount, MPIStructDef::getMPIType(&sendIDs[0]),
	       &recvIDs[0], &recvcounts[0], &displs[0], 
	       MPIStructDef::getMPIType(&sendIDs[0]), root, group.comm));
  
  if (grank == root) {
    // fill in the local array with all gathered data, after reordering them
    cf_assert(dtt->arraySize == totRecvcount);
    CFreal *const sarray = dtt->array;
    for (CFuint is = 0; is < totRecvcount; ++is) {
      cf_assert(is < recvIDs.size());
      cf_assert(is < recvbuf.size());
      sarray[recvIDs[is]] = recvbuf[is];  
    }
  } 
  
  CFLog(INFO, "StdConcurrentDataTransfer::gatherData() from namespace[" << nspSend 
	<< "] to namespace [" << nspRecv << "] => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void StdConcurrentDataTransfer::scatterData(const CFuint idx)
{ 
  // AL: have to involve only the ranks involved in this scattering
  SafePtr<DataToTrasfer> dtt = _socketName2data.find(_socketsSendRecv[idx]); 
  
  const string nspSend = dtt->nspSend;
  const string nspRecv = dtt->nspRecv;
  const string sendSocketStr = dtt->sendSocketStr;
  const string recvSocketStr = dtt->recvSocketStr;
  const string nspCoupling = getMethodData().getNamespace();
  
  CFLog(INFO, "StdConcurrentDataTransfer::scatterData() from namespace[" << nspSend 
	<< "] to namespace [" << nspRecv << "] => start\n");
  
  Group& group = PE::GetPE().getGroup(nspCoupling);
  const int rank  = PE::GetPE().GetRank("Default"); // rank in MPI_COMM_WORLD
  const int grank = PE::GetPE().GetRank(nspCoupling); // rank in coupling group
  const CFuint nbRanks = group.globalRanks.size();
    
  // build mapping from global to local DOF IDs on the receiving side
  if (PE::GetPE().isRankInGroup(rank, nspRecv)) {  
    if (_global2localIDs.size() == 0) {
      SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspRecv);
      if (_socketsConnType[idx] == "State") {
	fillMapGlobalToLocal<State*>(dtt, ds, _global2localIDs);
      }
      
      if (_socketsConnType[idx] == "Node") {
	fillMapGlobalToLocal<Node*>(dtt, ds, _global2localIDs);
      }
    }
  }
  
  vector<CFreal> sendbuf;
  vector<CFuint> sendIDs;
  vector<int> sendcounts(nbRanks, 0);
  vector<int> sendIDcounts(nbRanks, 0);
  
  // this scatters contributions from one rank in the "send" namespace 
  // to all ranks belonging to the "recv" namespace
  if (PE::GetPE().isRankInGroup(rank, nspSend)) {  
    SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspSend);
    if (_socketsConnType[idx] == "State") {
      fillSendCountsScatter<State*>(dtt, ds, sendcounts, sendIDcounts);
    }
    if (_socketsConnType[idx] == "Node")  {
      fillSendCountsScatter<Node*>(dtt, ds, sendcounts, sendIDcounts);
    }
  }
  
  int root = getRootProcess(nspSend, nspCoupling);
  
  MPIStruct msSizes;
  int ln[2];
  ln[0] = ln[1] = nbRanks;
  MPIStructDef::buildMPIStruct(&sendcounts[0], &sendIDcounts[0], ln, msSizes);
  
  MPIError::getInstance().check
    ("MPI_Bcast", "StdConcurrentDataTransfer::scatterData()", 
     MPI_Bcast(msSizes.start, 1, msSizes.type, root, group.comm));
    
  sendbuf.resize(sendcounts[nbRanks-1]);
  sendIDs.resize(sendIDcounts[nbRanks-1]);
  cf_assert(sendbuf.size() > 0);
  cf_assert(sendIDs.size() > 0);
  
  CFuint counter = 0;
  CFuint countID = 0;
  for (CFuint r = 0; r < nbRanks; ++r) {
    const CFuint sendSize = sendcounts[r]; 
    const CFuint sendIDSize = sendIDcounts[r]; 
    const CFuint stride = sendSize/sendIDSize;
    cf_assert(stride >= 1);
    
    if (grank == root) {
      CFreal *const dataToSend = dtt->array;
      cf_assert(dataToSend != CFNULL);
      for (CFuint s = 0; s < sendSize; ++s, ++counter) {
	sendbuf[s] = dataToSend[counter];
      }
      
      SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspSend);
      if (_socketsConnType[idx] == "State") {
	const string statesStr = nspSend + "_states";
	DataHandle<State*, GLOBAL> states = ds->getGlobalData<State*>(statesStr);
	for (CFuint s = 0; s < sendIDSize; ++s, ++countID) {
	  sendIDs[s] = states[countID]->getGlobalID(); 
	}
      }
      if (_socketsConnType[idx] == "Node")  {
	const string nodesStr = nspSend + "_nodes";
	DataHandle<Node*, GLOBAL> nodes = ds->getGlobalData<Node*>(nodesStr);
	for (CFuint s = 0; s < sendIDSize; ++s, ++countID) {
	  sendIDs[s] = nodes[countID]->getGlobalID(); 
	}
      }
    }
        
    MPIStruct ms;
    int ln[2];
    ln[0] = sendcounts[r];
    ln[1] = sendIDcounts[r];
    MPIStructDef::buildMPIStruct<CFreal, CFuint>(&sendbuf[0], &sendIDs[0], ln, ms);
    
    MPIError::getInstance().check
      ("MPI_Bcast", "StdConcurrentDataTransfer::scatterData()", 
       MPI_Bcast(ms.start, 1, ms.type, root, group.comm));
    
    if (grank != root) { 
      cf_assert(idx < _sendToRecvVecTrans.size());
      SafePtr<VarSetTransformer> sendToRecvTrans = _sendToRecvVecTrans[idx].getPtr();
      cf_assert(sendToRecvTrans.isNotNull());
      const CFuint sendStride = dtt->sendStride;
      const CFuint recvStride = dtt->recvStride;
      //  cf_assert(recvStride >= stride);
      RealVector tState(recvStride, static_cast<CFreal*>(NULL));
      RealVector state(sendStride, static_cast<CFreal*>(NULL));
      //  when current rank finds a globalID, it copies the data in corresponding localID position
      CFreal *const dataToRecv = dtt->array;
      cf_assert(dataToRecv != CFNULL);
      for (CFuint id = 0; id < sendIDSize; ++id) {
	bool found = false;
	const CFuint localID = _global2localIDs.find(sendIDs[id], found);
	if (found) {
	  // here you need a transformer
	  const CFuint startR = localID*recvStride;
	  const CFuint startS = id*stride;
	  cf_assert(stride == dtt->sendStride);
	  // state.wrap(sendStride, &dataToRecv[startR]); // recheck send/recv sizes here (startS instead?)
	  // tState.wrap(recvStride, &sendbuf[startS]);   // recheck send/recv sizes here (startR instead?)
	  state.wrap(sendStride, &sendbuf[startS]); 
	  tState.wrap(recvStride, &dataToRecv[startR]);
	  sendToRecvTrans->transform((const RealVector&)state, (RealVector&)tState);
	}
      }
    }
  }
  
  CFLog(INFO, "StdConcurrentDataTransfer::scatterData() from namespace[" << nspSend 
	<< "] to namespace [" << nspRecv << "] => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

int StdConcurrentDataTransfer::getRootProcess(const std::string& nsp, 
					      const std::string& nspCoupling) const
{
  const int rank  = PE::GetPE().GetRank("Default"); // rank in MPI_COMM_WORLD
  Group& group = PE::GetPE().getGroup(nspCoupling);
  
  int root = -1;
  int sendroot = -1; 
  if (PE::GetPE().isRankInGroup(rank, nsp)) {
    sendroot = PE::GetPE().GetRank(nspCoupling);
  }
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "StdConcurrentDataTransfer::getRootProcess()", 
     MPI_Allreduce(&sendroot, &root, 1, MPIStructDef::getMPIType(&root), 
		   MPI_MAX, group.comm));
  cf_assert(root >= 0);
  return root;
}
      
//////////////////////////////////////////////////////////////////////////////
   
void StdConcurrentDataTransfer::addDataToTransfer(const CFuint idx)
{
  DataToTrasfer* data = new DataToTrasfer();
  
  vector<string> sendRecv = StringOps::getWords(_socketsSendRecv[idx],'>');
  cf_assert(sendRecv.size() == 2);
  
  // namespace_socket (send)
  const string sendSocketStr = sendRecv[0];
  vector<string> nspSocketSend = StringOps::getWords(sendSocketStr,'_');
  cf_assert(nspSocketSend.size() == 2);
  
  // namespace_socket (recv)
  const string recvSocketStr = sendRecv[1];
  vector<string> nspSocketRecv = StringOps::getWords(recvSocketStr,'_');
  cf_assert(nspSocketRecv.size() == 2);
  
  const string nspSend    = nspSocketSend[0];
  const string socketSend = nspSocketSend[1];
  const string nspRecv    = nspSocketRecv[0];
  const string socketRecv = nspSocketRecv[1];
  
  CFLog(VERBOSE, "StdConcurrentDataTransfer::addDataToTransfer() => send: " 
	<< nspSend << "-" << socketSend << "\n");
  CFLog(VERBOSE, "StdConcurrentDataTransfer::addDataToTransfer() => recv: " 
	<< nspRecv << "-" << socketRecv << "\n");
  
  const int rank  = PE::GetPE().GetRank("Default"); // rank in MPI_COMM_WORLD
  Group& groupSend = PE::GetPE().getGroup(nspSend);
  Group& groupRecv = PE::GetPE().getGroup(nspRecv);
  const string sendLocal  = sendSocketStr + "_local";
  const string sendGlobal = sendSocketStr + "_global";
  const string recvLocal  = recvSocketStr + "_local";
  const string recvGlobal = recvSocketStr + "_global";
    
  data->nspSend = nspSend;
  data->nspRecv = nspRecv;
  data->sendSocketStr = sendSocketStr;
  data->recvSocketStr = recvSocketStr;
  data->nbRanksSend = groupSend.globalRanks.size();
  data->nbRanksRecv = groupRecv.globalRanks.size();
  
  vector<CFuint> sendRecvStridesIn(2, 0);
  
  // send data
  if (PE::GetPE().isRankInGroup(rank, nspSend)) {  
    SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspSend);
    
    // local data (CFreal)
    if (ds->checkData(sendSocketStr)) {
      DataHandle<CFreal> array = ds->getData<CFreal>(sendSocketStr);
      CFLog(VERBOSE, "P" << rank << " has socket " << sendSocketStr << "\n"); 
      
      CFuint dofsSize = 0;
      if (_socketsConnType[idx] == "State") {
	data->dofsName = nspSend + "_states";
	Framework::DataHandle<State*, GLOBAL> dofs = ds->getGlobalData<State*>(data->dofsName);
	dofsSize = dofs.size();
      }
      if (_socketsConnType[idx] == "Node") {
	data->dofsName = nspSend + "_nodes";
	Framework::DataHandle<Node*, GLOBAL> dofs = ds->getGlobalData<Node*>(data->dofsName);
	dofsSize = dofs.size();
      }
      
      data->array = &array[0]; 
      data->arraySize = array.size();
      cf_assert(data->arraySize > 0);
      sendRecvStridesIn[0] = array.size()/dofsSize;
    }
    // global data (State*)
    else if (ds->checkData(sendLocal) && ds->checkData(sendGlobal)) {
      CFLog(VERBOSE, "P" << rank << " has socket <State*> " << sendSocketStr << "\n"); 
      DataHandle<State*, GLOBAL> array = ds->getGlobalData<State*>(sendSocketStr);
      data->dofsName = sendSocketStr;
      data->array = array.getGlobalArray()->ptr();
      data->arraySize = array.size()*array[0]->size();
      cf_assert(data->arraySize > 0);
      sendRecvStridesIn[0] = array[0]->size();
    }
  }
  
  // recv data
  if (PE::GetPE().isRankInGroup(rank, nspRecv)) {  
    SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspRecv);
    
    // local data (CFreal)
    if (ds->checkData(recvSocketStr)) {
      DataHandle<CFreal> array = ds->getData<CFreal>(recvSocketStr);
      CFLog(VERBOSE, "P" << rank << " has socket " << recvSocketStr << "\n"); 
      
      CFuint dofsSize = 0;
      if (_socketsConnType[idx] == "State") {
	data->dofsName = nspRecv + "_states";
	Framework::DataHandle<State*, GLOBAL> dofs = ds->getGlobalData<State*>(data->dofsName);
	dofsSize = dofs.size();
      }
      if (_socketsConnType[idx] == "Node") {
	data->dofsName = nspRecv + "_nodes";
	Framework::DataHandle<Node*, GLOBAL> dofs = ds->getGlobalData<Node*>(data->dofsName);
	dofsSize = dofs.size();
      }
      
      data->array  = &array[0]; 
      data->arraySize = array.size();
      cf_assert(data->arraySize > 0);
      sendRecvStridesIn[1] = array.size()/dofsSize;
    }
    // global data (State*)
    else if (ds->checkData(recvLocal) && ds->checkData(recvGlobal)) {
      CFLog(VERBOSE, "P" << rank << " has socket <State*> " << recvSocketStr << "\n"); 
      DataHandle<State*, GLOBAL> array = ds->getGlobalData<State*>(recvSocketStr);
      CFLog(VERBOSE, "P" << rank << " has socket " << recvSocketStr << " with sizes = [" 
	    << array.getLocalSize() << ", " << array.getGlobalSize() << "]\n"); 
      cf_assert(array.getLocalSize() == array.getGlobalSize());
      cf_assert(array.size() == array.getGlobalSize());
      
      data->dofsName = recvSocketStr;
      data->array = array.getGlobalArray()->ptr();
      data->arraySize = array.size()*array[0]->size();
      cf_assert(data->arraySize > 0);
      sendRecvStridesIn[1] = array[0]->size();
    }
  }
  
  vector<CFuint> sendRecvStridesOut(2,0);
  const string nspCoupling = getMethodData().getNamespace();
  Group& group = PE::GetPE().getGroup(nspCoupling);
  
  MPIError::getInstance().check
  ("MPI_Allreduce", "StdConcurrentDataTransfer::addDataToTransfer()", 
   MPI_Allreduce(&sendRecvStridesIn[0], &sendRecvStridesOut[0], 2,
		 MPIStructDef::getMPIType(&sendRecvStridesIn[0]), MPI_MAX, group.comm));
  
  data->sendStride = sendRecvStridesOut[0];
  data->recvStride = sendRecvStridesOut[1];
  cf_assert(data->sendStride > 0);
  cf_assert(data->recvStride > 0);
  // cf_assert(data->sendStride >= data->recvStride);
  
  _socketName2data.insert(_socketsSendRecv[idx], data);
}
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD
