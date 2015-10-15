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
#include "ConcurrentCoupler/StdConcurrentReduce.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdConcurrentReduce, 
		      ConcurrentCouplerData, 
		      ConcurrentCouplerModule> 
stdConcurrentReduceProvider("StdConcurrentReduce");

//////////////////////////////////////////////////////////////////////////////

void StdConcurrentReduce::defineConfigOptions(Config::OptionList& options)
{  
  options.addConfigOption< vector<string> >
    ("SocketsSendRecv","Sockets to transfer, for example: Namespace1_send>Namespace2_recv (no space on both sides of \">\".");
  options.addConfigOption< vector<string> >
    ("SocketsConnType","Connectivity type for sockets to transfer (State or Node): this is ne1eded to define global IDs.");
  options.addConfigOption< vector<string> >
    ("Operation","Name of operation to perform during the reduction for each socket (SUM, SUB, MIN. MAX, PROD).");
}
      
//////////////////////////////////////////////////////////////////////////////

StdConcurrentReduce::StdConcurrentReduce(const std::string& name) :
  ConcurrentCouplerCom(name),
  m_sockets(),
  socket_states("states"),
  m_socketName2data(),
  m_nameToMPIOp(),
  m_isReduceRank(),
  m_reduceNsp("")
{
  addConfigOptionsTo(this);
  
  m_nameToMPIOp.insert("SUM",  MPI_SUM);
  m_nameToMPIOp.insert("PROD", MPI_PROD);
  m_nameToMPIOp.insert("MAX",  MPI_MAX);
  m_nameToMPIOp.insert("MIN",  MPI_MIN);
  
  m_socketsSendRecv = vector<string>();
  setParameter("SocketsSendRecv",&m_socketsSendRecv);
  
  m_socketsConnType = vector<string>();
  setParameter("SocketsConnType",&m_socketsConnType);
  
  m_operation = vector<string>();
  setParameter("Operation",&m_operation);
}
      
//////////////////////////////////////////////////////////////////////////////

StdConcurrentReduce::~StdConcurrentReduce()
{
  for (CFuint i =0 ; i < m_socketName2data.size(); ++i) {
    delete m_socketName2data[i];
  }  
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdConcurrentReduce::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = 
    m_sockets.getAllSinkSockets();
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdConcurrentReduce::configure ( Config::ConfigArgs& args )
{
  ConcurrentCouplerCom::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////
  
void StdConcurrentReduce::setup()
{
  // "SUM" is default operation if user does not specify the correct number of 
  // operations, once per socket
  if (m_operation.size() != m_socketsSendRecv.size()) {
    m_operation.resize(m_socketsSendRecv.size(), "SUM");
  }
}
      
//////////////////////////////////////////////////////////////////////////////
      
void StdConcurrentReduce::execute()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "StdConcurrentReduce::execute() => start\n");
  
  // this should go in the setup, but it uses blocking MPI collective calls 
  // here is less harmful 
  if (m_socketName2data.size() == 0) {
    createReduceGroup();
    
    if (isReductionRank()) {
      for (CFuint i = 0; i < m_socketsSendRecv.size(); ++i) {
	addDataToTransfer(i);
      }
    }
  } 
  
  if (isReductionRank()) {
    // Commands "live" in their Method's namespace, therefore only ranks belonging to 
    // nspCoupling (in this case) get here, because all other ranks are filtered out
    // at the higher level (in the corresponding Method's member function)
    Group& group = PE::GetPE().getGroup(getName());
    
    vector<CFreal> recvBuf;
    for (CFuint i = 0; i < m_socketsSendRecv.size(); ++i) {
      SafePtr<DataToTrasfer> dtt = m_socketName2data.find(m_socketsSendRecv[i]); 
      // note this limits the maximum size to be transferred to an array of size ~2 billion 
      const CFuint arraySize = dtt->arraySize;
      recvBuf.resize(arraySize, 0.);
      
      CFreal* sendBuf = dtt->array;
      cf_assert(sendBuf != CFNULL);
      
      MPIError::getInstance().check
	("MPI_Allreduce", "StdConcurrentReduce::execute()", 
	 MPI_Allreduce(&sendBuf[0], &recvBuf[0], (int)arraySize, 
		       MPIStructDef::getMPIType(&recvBuf[0]), dtt->operation, group.comm));
      
      // overwrite the socket with the total value in each processor
      for (CFuint s = 0; s < arraySize; ++s) {
	sendBuf[s] = recvBuf[s];
      }
    }
  }
  
  CFLog(INFO, "StdConcurrentReduce::execute() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////
   
 void StdConcurrentReduce::addDataToTransfer(const CFuint idx)
 {
   CFLog(VERBOSE, "StdConcurrentReduce::addDataToTransfer() for [" << idx << "] => start\n");
   
   const int rank  = PE::GetPE().GetRank("Default");
   DataToTrasfer* data = new DataToTrasfer();
   CFuint sendRecvStridesIn = 0;
   
   // reduction is only applied to ranks belonging to the group whose name is getName()
   const size_t found1 = m_socketsSendRecv[idx].find(">");
   if (found1 == string::npos) {
     // if the ">" is not found, separate namespace from socket name
     vector<string> nsp2Socket = StringOps::getWords(m_socketsSendRecv[idx], '_');
     const string nspSend      = m_reduceNsp;
     cf_assert(nspSend != "");
     const string sendSocketStr = nspSend + "_" + nsp2Socket[1];
     const string sendLocal  = sendSocketStr + "_local";
     const string sendGlobal = sendSocketStr + "_global";
     SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspSend);
     data->nspSend = nspSend;
     data->sendSocketStr = sendSocketStr;
     // data->nbRanksSend = groupSend.globalRanks.size();
     
     // local data (CFreal)
     if (ds->checkData(sendSocketStr)) {
       DataHandle<CFreal> array = ds->getData<CFreal>(sendSocketStr);
       CFLog(VERBOSE, "P" << rank << " has socket " << sendSocketStr << "\n"); 
       
       CFuint dofsSize = 0;
       if (m_socketsConnType[idx] == "State") {
	 data->dofsName = nspSend + "_states";
	 DataHandle<State*, GLOBAL> dofs = ds->getGlobalData<State*>(data->dofsName);
	 dofsSize = dofs.size();
       }
       if (m_socketsConnType[idx] == "Node") {
	 data->dofsName = nspSend + "_nodes";
	 DataHandle<Node*, GLOBAL> dofs = ds->getGlobalData<Node*>(data->dofsName);
	 dofsSize = dofs.size();
       }
       
       data->array = &array[0]; 
       data->arraySize = array.size();
       cf_assert(data->arraySize > 0);
       sendRecvStridesIn = array.size()/dofsSize;
     }
     // global data (State*)
     else if (ds->checkData(sendLocal) && ds->checkData(sendGlobal)) {
       CFLog(VERBOSE, "P" << rank << " has socket <State*> " << sendSocketStr << "\n"); 
       DataHandle<State*, GLOBAL> array = ds->getGlobalData<State*>(sendSocketStr);
       data->dofsName = sendSocketStr;
       data->array = array.getGlobalArray()->ptr();
       data->arraySize = array.size()*array[0]->size();
       cf_assert(data->arraySize > 0);
       sendRecvStridesIn = array[0]->size();
     }
     
     CFuint sendRecvStridesOut = 0;
     Group& group = PE::GetPE().getGroup(getName());
     
     cf_assert(isReductionRank());
     
     MPIError::getInstance().check
       ("MPI_Allreduce", "StdConcurrentDataTransfer::addDataToTransfer()", 
	MPI_Allreduce(&sendRecvStridesIn, &sendRecvStridesOut, 1,
		      MPIStructDef::getMPIType(&sendRecvStridesIn), MPI_MAX, group.comm));
     
     cf_assert(sendRecvStridesIn == sendRecvStridesOut);
     
     data->sendStride = sendRecvStridesOut; 
     CFLog(VERBOSE, "P" << rank << " has data->sendStride = " << data->sendStride << "\n"); 
     cf_assert(data->sendStride > 0);
     cf_assert(idx < m_operation.size());
     data->operation = m_nameToMPIOp.find(m_operation[idx]); 
     m_socketName2data.insert(m_socketsSendRecv[idx], data);
   }
   else {
     throw NotImplementedException
       (FromHere(), "StdConcurrentReduce::addDataToTransfer() => socket name has \">\"");
   }
   
   CFLog(VERBOSE, "StdConcurrentReduce::addDataToTransfer() for [" << idx << "] => end\n");
 }
      
//////////////////////////////////////////////////////////////////////////////
 
void StdConcurrentReduce::createReduceGroup()
{
  CFLog(VERBOSE, "StdConcurrentReduce::createReduceGroup() => start\n");
  
  const string nspCoupling = getMethodData().getNamespace();
  const Group& nspGroup = PE::GetPE().getGroup(nspCoupling);
  const CFuint nspRanksSize = nspGroup.globalRanks.size();
  const int nspRank  = PE::GetPE().GetRank(nspCoupling);
  cf_assert(nspRank < nspRanksSize);
  
  vector<int> isReduceRank(nspRanksSize, 0); 
  m_isReduceRank.resize(nspRanksSize, 0);
  
  // AL: we suppose that each socket in m_socketsSendRecv belongs to the same 
  // namespace, therefore we only test the existence of the first socket name  
  const string testSocket = m_socketsSendRecv[0]; 
  // the socket name must not include ">"
  const size_t found1 = testSocket.find(">"); 
  if (found1 == string::npos) {
    // if the ">" is not found, separate namespace from socket name
    vector<string> nsp2Socket = StringOps::getWords(testSocket, '_');
    const string nsp          = nsp2Socket[0];
    const string socketSend   = nsp2Socket[1];
    
    // scan all the ranks to check if they contain a socket called
    // nsp + rank + "_" + socketSend  (e.g. Namespace0_data, Namespace1_data, etc.)
    for (CFuint rk = 0; rk < nspRanksSize; ++rk) {
      const string nspSend = nsp + StringOps::to_str(rk);
      if (PE::GetPE().checkGroup(nspSend)) {
	const string sendSocket = nspSend + "_" + socketSend;
	const string sendLocal  = sendSocket + "_local";
	const string sendGlobal = sendSocket + "_global";
	SafePtr<DataStorage> ds = getMethodData().getDataStorage(nspSend);
	
	// local data (CFreal)
	if ((ds->checkData(sendSocket) && 
	     ds->getData<CFreal>(sendSocket).size()>0) ||
	    
	    (ds->checkData(sendLocal) && ds->checkData(sendGlobal) && 
	     ds->getGlobalData<State*>(sendSocket).size()>0)) {
	  isReduceRank[nspRank] = 1;
	  m_reduceNsp = nspSend;
	  break;
	}
      }
    }
    
    MPIError::getInstance().check
      ("MPI_Allreduce", "StdConcurrentReduce::createReduceGroup()", 
       MPI_Allreduce(&isReduceRank[0], &m_isReduceRank[0], nspRanksSize, 
		     MPIStructDef::getMPIType(&isReduceRank[0]), MPI_MAX, nspGroup.comm));
    
    CFLog(VERBOSE, "StdConcurrentReduce::createReduceGroup() => P[" << PE::GetPE().GetRank("Default") 
	  << "] in \"Default\", P[" << nspRank << "] in " << "\"" << nspCoupling << "\" in(1) or out(0) : " 
	  << m_isReduceRank[nspRank] << " in \"" << m_reduceNsp << "\"\n");
    
    vector<int> ranks;
    for (int rk = 0; rk < nspRanksSize; ++rk) {
      if (m_isReduceRank[rk] == 1) {ranks.push_back(rk);}
    }
    cf_assert(ranks.size() > 0);
    
    // here we create a subgroup of the current coupling namespace 
    PE::GetPE().createGroup(nspCoupling, getName(), ranks, true);
    
    const string msg = "StdConcurrentReduce::createReduceGroup() => Ranks for group [" + getName() + "] = ";
    CFLog(VERBOSE, CFPrintContainer<vector<int> >(msg, &ranks));
  }
  
  CFLog(VERBOSE, "StdConcurrentReduce::createReduceGroup() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////
 
    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD
