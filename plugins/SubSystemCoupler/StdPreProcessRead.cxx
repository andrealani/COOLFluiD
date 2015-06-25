#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/VolumeIntegratorImpl.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/StdPreProcessRead.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdPreProcessRead, SubSysCouplerData, SubSystemCouplerModule> stdPreProcessReadProvider("StdPreProcessRead");

//////////////////////////////////////////////////////////////////////////////

StdPreProcessRead::StdPreProcessRead(const std::string& name) : CouplerCom(name),
  _sockets(),
  _nbOtherStates(0),
  _alreadyReadSockets(0),
  socket_states("states"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdPreProcessRead::~StdPreProcessRead()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessRead::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

    CFLog(INFO, "Command [" << getName()
             << "] Data namespace [" << getMethodData().getNamespace()
             << "]\n");

  CouplerCom::configure(args);

  _sockets.setParentNamespace( getMethodData().getNamespace() );

    CFLog(INFO, "Command [" << getName()
             << "] Data namespace [" << getMethodData().getNamespace()
             << "]\n");

  typedef std::vector<Common::SafePtr<CommandGroup> > InterfaceList;
  InterfaceList interfaces = getMethodData().getInterfaces();

  try {
    const std::string nsp = getMethodData().getNamespace();
  InterfaceList::iterator itr = interfaces.begin();
  for(; itr != interfaces.end(); ++itr) {

    const std::vector<std::string>& comNames = (*itr)->getComNames();

CF_DEBUG_POINT;

    // check if this command applies to this Interface
    if(count(comNames.begin(),comNames.end(),getName()) != 0) {

      const std::string interfaceName = (*itr)->getName();

      for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++iProc)
      {
        const vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);
        for (CFuint iTRS = 0; iTRS < otherTrsNames.size(); ++iTRS)
        {
          /// Create the datahandle for the otherSubSystem Coordinates
          const std::string otherTrsName = otherTrsNames[iTRS];

          /// Get the datahandle for the coordinates of the OTHER SubSystem
          const vector<std::string> socketCoordNames =
            getMethodData().getOtherCoupledCoordName(interfaceName,otherTrsName,iProc);
          const vector<std::string> socketAcceptNames =
            getMethodData().getOtherCoupledAcceptedName(interfaceName,otherTrsName,iProc);
          const vector<std::string> socketDataNames =
            getMethodData().getOtherCoupledDataName(interfaceName,otherTrsName,iProc);
CF_DEBUG_POINT;
          for(CFuint iType=0;iType < socketCoordNames.size();iType++)
          {
            _sockets.createSocketSource<RealVector>(socketCoordNames[iType]);
            _sockets.createSocketSource<CFreal>(socketAcceptNames[iType]);
            _sockets.createSocketSource<RealVector>(socketDataNames[iType]);
CF_DEBUG_POINT;
          }
        } //loop over processors
      } // other trs
    } // if check
  } // interfaces

  } catch (Exception& e)
  {
    CFout << e.what() << "\n" << CFendl;
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessRead::executeOnTrs()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "Reading Coordinates of coupled coordinates for TRS" <<getCurrentTRS()->getName() << "\n");
  
  const std::string nsp = getMethodData().getNamespace();
  const bool isParallel = Common::PE::GetPE().IsParallel();
  if(isParallel)
  {
    Common::PE::GetPE().setBarrier(nsp);
    
    for (CFuint i = 0; i < Common::PE::GetPE().GetProcessorCount(nsp); ++i)
    {
      if (i == Common::PE::GetPE().GetRank(nsp))
      {
        for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++iProc)
        {
          executeReadOnTrs(iProc);
        }
      }
      Common::PE::GetPE().setBarrier(nsp);
    }
  }
  else
  {
    executeReadOnTrs(0);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessRead::executeReadOnTrs(const CFuint iProc)
{
  CFAUTOTRACE;

  const std::string interfaceName = getCommandGroupName();
  const vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);

  for (CFuint iTRS = 0; iTRS< otherTrsNames.size(); ++iTRS)
  {
    const vector<std::string> socketCoordNames =
      getMethodData().getOtherCoupledCoordName(interfaceName,otherTrsNames[iTRS],iProc);
    const vector<std::string> socketAcceptNames =
      getMethodData().getOtherCoupledAcceptedName(interfaceName,otherTrsNames[iTRS],iProc);
    const vector<std::string> socketDataNames =
      getMethodData().getOtherCoupledDataName(interfaceName,otherTrsNames[iTRS],iProc);

    for(CFuint iType=0; iType < socketCoordNames.size();iType++)
    {
      DataHandle< RealVector> interfaceCoords =
        _sockets.getSocketSource<RealVector>(socketCoordNames[iType])->getDataHandle();
      DataHandle< CFreal> isAccepted =
        _sockets.getSocketSource<CFreal>(socketAcceptNames[iType])->getDataHandle();
      DataHandle< RealVector> interfaceData =
        _sockets.getSocketSource<RealVector>(socketDataNames[iType])->getDataHandle();

      const std::string socketName = socketCoordNames[iType];

      if (std::count(_alreadyReadSockets.begin(),_alreadyReadSockets.end(),socketName) == 0) {
        // Temporary: transfer info through files...
        // Read the otherSubSystem_COORD files
        readFile(socketName);

        _alreadyReadSockets.push_back(socketName);
      }

      /// resize the datahandles
      const CFuint nbCoupledStates = interfaceCoords.size();
      isAccepted.resize(nbCoupledStates);
      interfaceData.resize(nbCoupledStates);

      ///Resize the data to the correct size
      const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
      const CFuint transferedSize = nbEqs;
      const CFuint transformedSize =
        getMethodData().getPreVariableTransformer(interfaceName)->getTransformedSize(transferedSize);

      for(CFuint iData=0;iData< interfaceData.size(); iData++)
      {
        interfaceData[iData].resize(transformedSize);
      }
    }
  } // other trs name
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessRead::readFile(const std::string& socketName)
{
  CFAUTOTRACE;

  DataHandle<RealVector> interfaceCoord =
    _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

  boost::filesystem::path fname =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(socketName);

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& finCoord = fhandle->open(fname);

  CFuint lineNb = 0;
  CFuint nbOtherStates;
  std::string line;
  vector<std::string> words;

  // read nb states
  getWordsFromLine(finCoord,line,lineNb,words);
  cf_assert(words.size() == 1);
  nbOtherStates = Common::StringOps::from_str<CFint>(words[0]);

  // Resize DataHandle
  interfaceCoord.resize(nbOtherStates);
  for (CFuint i=0; i< nbOtherStates; ++i)
  {
  getWordsFromLine(finCoord,line,lineNb,words);
  // Resize to the dimension of the problem
  (interfaceCoord[i]).resize(words.size());
//  CFout << "Size of variables: " << words.size() << "\n";
    for (CFuint j=0; j<words.size();++j){
      (interfaceCoord[i])[j] =  Common::StringOps::from_str<CFreal>(words[j]);
    }
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessRead::getWordsFromLine(ifstream& fin,
                                     std::string& line,
                                     CFuint&  lineNb,
                                     vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdPreProcessRead::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = _sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdPreProcessRead::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();

  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD
