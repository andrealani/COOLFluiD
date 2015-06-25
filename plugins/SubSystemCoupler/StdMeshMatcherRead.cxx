

#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "SubSystemCoupler/StdMeshMatcherRead.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"

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

MethodCommandProvider<StdMeshMatcherRead, SubSysCouplerData, SubSystemCouplerModule> stdMeshMatcherReadProvider("StdMeshMatcherRead");

//////////////////////////////////////////////////////////////////////////////

StdMeshMatcherRead::StdMeshMatcherRead(const std::string& name) :
  CouplerCom(name),
  _sockets()
{
}

//////////////////////////////////////////////////////////////////////////////

StdMeshMatcherRead::~StdMeshMatcherRead()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdMeshMatcherRead::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdMeshMatcherRead::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = _sockets.getAllSourceSockets();;

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherRead::configure ( Config::ConfigArgs& args )
{
  CouplerCom::configure(args);

  typedef std::vector<Common::SafePtr<CommandGroup> > InterfaceList;
  InterfaceList interfaces = getMethodData().getInterfaces();

  try {

  InterfaceList::iterator itr = interfaces.begin();
  for(; itr != interfaces.end(); ++itr) {

    const std::vector<std::string>& comNames = (*itr)->getComNames();

    // check if this command applies to this Interface
    if(count(comNames.begin(),comNames.end(),getName()) != 0) {

      //Read Part
      const std::string interfaceName = (*itr)->getName();
      const std::vector<std::string> trsNames = (*itr)->getTrsNames();

      std::vector<std::string>::const_iterator trs = trsNames.begin();
      for(; trs != trsNames.end(); ++trs)
      {
        const std::string trsName = *trs;

        /// Get the datahandle for the acceptance flags to be transfered TO THE CURRENT SubSystem
        const vector<std::string> socketAcceptNames =
          getMethodData().getThisCoupledAcceptedName(interfaceName,trsName);
        const vector<std::string> socketDataNames =
          getMethodData().getThisCoupledDataName(interfaceName,trsName);
        for(CFuint iType=0;iType < socketDataNames.size();iType++)
        {
          _sockets.createSocketSink<CFreal>(socketAcceptNames[iType]);
          _sockets.createSocketSink<RealVector>(socketDataNames[iType]);
          _sockets.createSocketSink<RealVector>(socketDataNames[iType] + "_PAST");
          _sockets.createSocketSink<RealVector>(socketDataNames[iType] + "_ORIGINAL");

          _sockets.createSocketSource<CFuint>(socketAcceptNames[iType] + "PAR");
        }
      } // trs
    } // if check
  } // interfaces

  } catch (Exception& e)
  {
    CFout << e.what() << "\n" << CFendl;
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherRead::execute()
{

  CFAUTOTRACE;
  const std::string nsp = getMethodData().getNamespace();
  //Here no need for barrier because each processor reads different files
  for (_iProc = 0; _iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++_iProc) {
    executeRead();
  }

  resizeDataHandles();
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherRead::executeRead()
{
  CFAUTOTRACE;

  const std::string interfaceName = getCommandGroupName();

  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    const std::string currentTrsName = getTrsName(iTRS);

    ///Read the file and put it into the datahandle
    const vector<std::string> localAcceptedSocketNames =
      getMethodData().getThisCoupledAcceptedName(interfaceName, currentTrsName);

    ///File name
    const vector<std::string> parAcceptedSocketNames =
      getMethodData().getThisCoupledAcceptedName(interfaceName, currentTrsName, _iProc);

    for(CFuint iType=0;iType< localAcceptedSocketNames.size();iType++)
    {
      readIsAcceptedFile(parAcceptedSocketNames[iType]);
      assembleIsAcceptedFiles(localAcceptedSocketNames[iType]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherRead::resizeDataHandles()
{
  const std::string interfaceName = getCommandGroupName();

  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    const std::string currentTrsName = getTrsName(iTRS);

    const vector<std::string> localAcceptedSocketNames =
      getMethodData().getThisCoupledAcceptedName(interfaceName, currentTrsName);
    const vector<std::string> dataSocketNames =
      getMethodData().getThisCoupledDataName(interfaceName,currentTrsName);

    cf_assert(dataSocketNames.size() == localAcceptedSocketNames.size());
    for(CFuint iType=0;iType< localAcceptedSocketNames.size();iType++)
    {
      //first compute the total number of accepted states

      DataHandle< CFreal> isAccepted =
        _sockets.getSocketSink<CFreal>(localAcceptedSocketNames[iType])->getDataHandle();
      CFuint nbAcceptedStates = 0;
      for(CFuint iState=0; iState < isAccepted.size() ;++iState){
        if(isAccepted[iState] >= 0.)
        {
          CFLog (VERBOSE, "State["<<iState <<"] accepted with distance: "<< isAccepted[iState]<< "\n" );
          ++nbAcceptedStates;
        }
      }
CFout << "MeshMatching nbAcceptedStates/nbLocalStates: "<<nbAcceptedStates <<"/" <<isAccepted.size() << "\n";

      //Resize the Parallel Data Index and fill it
      DataHandle< CFuint> isAcceptedParallelIndex =
        _sockets.getSocketSource<CFuint>(localAcceptedSocketNames[iType] + "PAR")->getDataHandle();
      isAcceptedParallelIndex.resize(nbAcceptedStates);

      CFuint idx=0;
      for(CFuint iState=0; iState < isAccepted.size() ;++iState){
        if(isAccepted[iState] >= 0.){
          isAcceptedParallelIndex[idx] = _tempParallelIndex[iState];
          ++idx;
        }
      }

      //Resize the Data datahandles
      const std::string dataSocketName = dataSocketNames[iType];
      DataHandle< RealVector> interfaceData =
        _sockets.getSocketSink<RealVector>(dataSocketName)->getDataHandle();
      DataHandle< RealVector> interfacePastData =
        _sockets.getSocketSink<RealVector>(dataSocketName + "_PAST")->getDataHandle();
      DataHandle< RealVector> interfaceOrigData =
        _sockets.getSocketSink<RealVector>(dataSocketName + "_ORIGINAL")->getDataHandle();

      interfaceData.resize(nbAcceptedStates);
      interfacePastData.resize(nbAcceptedStates);
      interfaceOrigData.resize(nbAcceptedStates);

      const CFuint transferedSize = getMethodData().getThisTransferedSize(interfaceName);
      const CFuint transformedSize =
        getMethodData().getPostVariableTransformer(interfaceName)->getTransformedSize(transferedSize);

      for(CFuint iData=0;iData< interfaceData.size(); iData++)
      {
        interfaceOrigData[iData].resize(transferedSize);
        interfaceData[iData].resize(transformedSize);
        interfacePastData[iData].resize(transformedSize);
      }
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherRead::readIsAcceptedFile(const std::string dataHandleName)
{
  CFAUTOTRACE;

  boost::filesystem::path fname =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(dataHandleName);

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(fname);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;

  // read nb states
  getWordsFromLine(fin,line,lineNb,words);
  cf_assert(words.size() == 1);
  const CFuint nbStates = Common::StringOps::from_str<CFint>(words[0]);

  // Resize DataHandle
  _tempIsAccepted.resize(nbStates);
  for (CFuint iState=0; iState < nbStates; ++iState)
  {
    getWordsFromLine(fin,line,lineNb,words);

    _tempIsAccepted[iState] = Common::StringOps::from_str<CFreal>(words[0]);
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherRead::assembleIsAcceptedFiles(const std::string localDataHandleName)
{
  CFAUTOTRACE;

  DataHandle< CFreal> localIsAccepted =
    _sockets.getSocketSink<CFreal>(localDataHandleName)->getDataHandle();

  //Resize datahandle the first time
  const CFuint nbStates = _tempIsAccepted.size();
  if(localIsAccepted.size() == 0){
    localIsAccepted.resize(nbStates);
    _tempParallelIndex.resize(nbStates);
    for (CFuint iState=0; iState < nbStates; ++iState){
      localIsAccepted[iState] = -1.;
      _tempParallelIndex[iState] = 0;
    }
  }

  for (CFuint iState=0; iState < nbStates; ++iState)
  {
    if((_tempIsAccepted[iState]>= 0.) &&
       ((_tempIsAccepted[iState] < localIsAccepted[iState]) || (localIsAccepted[iState] < 0.) ))
    {
      localIsAccepted[iState] = _tempIsAccepted[iState];
      _tempParallelIndex[iState] = _iProc;
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherRead::getWordsFromLine(ifstream& fin,
                                      std::string& line,
                                      CFuint&  lineNb,
                                      vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace MeshMatcher

  } // namespace Numerics

} // namespace COOLFluiD
