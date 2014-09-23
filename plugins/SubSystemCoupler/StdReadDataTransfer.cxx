
#include <boost/filesystem/path.hpp>

#include "Common/FilesystemException.hh"
#include "Common/NotImplementedException.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodCommand.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/VolumeIntegratorImpl.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/StdReadDataTransfer.hh"
#include "Framework/NamespaceSwitcher.hh"

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

MethodCommandProvider<StdReadDataTransfer, SubSysCouplerData, SubSystemCouplerModule> StdReadDataTransferProvider("StdReadDataTransfer");

//////////////////////////////////////////////////////////////////////////////

StdReadDataTransfer::StdReadDataTransfer(const std::string& name) :
  CouplerCom(name),
  _sockets(),
  socket_states("states"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdReadDataTransfer::~StdReadDataTransfer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdReadDataTransfer::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();;

  result.push_back(&socket_states);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::configure ( Config::ConfigArgs& args )
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

      const std::string interfaceName = (*itr)->getName();
      const std::vector<std::string> trsNames = (*itr)->getTrsNames();

      std::vector<std::string>::const_iterator trs = trsNames.begin();
      for(; trs != trsNames.end(); ++trs) {
        const std::string trsName = *trs;

        // place socket for accepted data
        const vector<std::string> socketAcceptNames =
          getMethodData().getThisCoupledAcceptedName(interfaceName,trsName);
        const vector<std::string> socketDataNames =
          getMethodData().getThisCoupledDataName(interfaceName,trsName);
        const vector<std::string> socketCoordNames =
          getMethodData().getThisCoupledCoordName(interfaceName,trsName);
        const vector<std::string> socketConnectNames =
          getMethodData().getThisCoupledConnectedName(interfaceName,trsName);
        for(CFuint iType=0;iType < socketAcceptNames.size();iType++)
        {
          _sockets.createSocketSink<CFreal>(socketAcceptNames[iType]);
          _sockets.createSocketSink<RealVector>(socketDataNames[iType]);
          _sockets.createSocketSink<RealVector>(socketDataNames[iType] + "_PAST");
          _sockets.createSocketSink<RealVector>(socketDataNames[iType] + "_ORIGINAL");
          _sockets.createSocketSink<RealVector>(socketCoordNames[iType]);
          _sockets.createSocketSink<std::vector<SubSysCouplerData::GeoEntityIdx> >(socketConnectNames[iType]);

          _sockets.createSocketSink<CFuint>(socketAcceptNames[iType] + "PAR");
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

void StdReadDataTransfer::setup()
{
  CouplerCom::setup();

  _interfaceName = getCommandGroupName();
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  ///@todo fix this for parallel
  if(!(Common::PE::GetPE().IsParallel()))
  {
    vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
    for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
    {
      _currentTrsName = getTrsName(iTRS);

      vector<std::string> socketDataNames =
        getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName);
      vector<std::string> coordTypes = getMethodData().getThisCoordType(_interfaceName);

      for(CFuint iType=0;iType < socketDataNames.size();iType++)
      {
        prepareNormFile(socketDataNames[iType]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::execute()
{
  CFAUTOTRACE;

  for (_iProc = 0; _iProc < Common::PE::GetPE().GetProcessorCount(); ++_iProc)
  {
    executeRead();
  }

  transformReceivedData();

}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::executeRead()
{

  CFAUTOTRACE;

  _interfaceName = getCommandGroupName();
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    _currentTrsName = getTrsName(iTRS);

    vector<std::string> socketAcceptedNames = getMethodData().getThisCoupledAcceptedName(_interfaceName,_currentTrsName);
    vector<std::string> socketDataNames = getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName);
    vector<std::string> parDataFileNames = getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName,_iProc);
    vector<std::string> parAcceptedFileNames = getMethodData().getThisCoupledAcceptedName(_interfaceName,_currentTrsName,_iProc);
    vector<std::string> coordTypes = getMethodData().getThisCoordType(_interfaceName);

    for(CFuint iType=0;iType < socketDataNames.size();iType++)
    {
      ///Read the file and put it into the datahandle
      if(getMethodData().isTransferFiles())
      {
        readFile(parDataFileNames[iType],parAcceptedFileNames[iType],socketDataNames[iType], socketAcceptedNames[iType]);
      }
      else{
        const bool isParallel = Common::PE::GetPE().IsParallel();
        cf_assert(!isParallel);
        readFromDataHandle(socketDataNames[iType]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::transformReceivedData()
{

  CFAUTOTRACE;

  _interfaceName = getCommandGroupName();
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  for (CFuint iTRS=0; iTRS < trs.size(); iTRS++)
  {
    _currentTrsName = getTrsName(iTRS);
    const vector<std::string> socketDataNames =
      getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName);
    const vector<std::string> coordTypes =
      getMethodData().getThisCoordType(_interfaceName);

    for(CFuint iType=0;iType < socketDataNames.size();iType++)
    {

      varTransfo->setCurrentInterface(_interfaceName, _currentTrsName, coordTypes[iType]);
      if(coordTypes[iType] == "Nodes")
      {
        transformReceivedNodalData();
      }
      if(coordTypes[iType] == "States")
      {
        transformReceivedStatesData();
      }
      if(coordTypes[iType] == "Gauss")
      {
        transformReceivedQuadratureData();
      }
      if(coordTypes[iType] == "Ghost")
      {
        transformReceivedGhostData();
      }

//      CFout << "Interface: " << _interfaceName <<"   -   TRS: "<< _currentTrsName << "\n";
      const bool isParallel = Common::PE::GetPE().IsParallel();
      ///@todo make it work in parallel
      if(!isParallel)
      {
        computeNorm(socketDataNames[iType]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::transformReceivedQuadratureData()
{
  ///Get Received DataHandle
  const std::string socketDataName = getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName,"Gauss");
  DataHandle< RealVector> interfaceData =
      _sockets.getSocketSink<RealVector>(socketDataName)->getDataHandle();
  DataHandle< RealVector> originalData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_ORIGINAL")->getDataHandle();

  DataHandle< RealVector> pastData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_PAST")->getDataHandle();

  const std::string socketAcceptedName = getMethodData().getThisCoupledAcceptedName(_interfaceName,_currentTrsName,"Gauss");
  DataHandle< CFreal> acceptedData =
      _sockets.getSocketSink< CFreal>(socketAcceptedName)->getDataHandle();

  /// Transform before setting the value to the DataHandle
  /// Call a variableTransformer
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  Common::SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_currentTrsName);
  geoData.trs = trs;
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  RealVector* tOtherState(CFNULL);
  RealVector coord(PhysicalModelStack::getActive()->getDim());

  ///Loop over geo entities to fill the datahandle
  //index going through the whole TRS
  CFuint idx = 0;
  // index going through the accepted points in the TRS
  CFuint acceptedIdx = 0;
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
  {
    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity* currFace = geoBuilder->buildGE();

    vector<SubSysCouplerData::GeoEntityIdx> faces(1);
    faces[0].first = trs;
    faces[0].second = geoData.idx;

    // gets the correct integrator corresponding to the
    // shape of the element where you want to integrate
    VolumeIntegratorImpl *const impl =
      getMethodData().getCollaborator<SpaceMethod>()->getVolumeIntegrator()->getSolutionIntegrator(currFace);

    std::vector<RealVector> quadraturePoints = impl->getQuadraturePointsCoordinates();
    for (CFuint iQuad = 0; iQuad < quadraturePoints.size(); ++iQuad)
    {
      coord = currFace->computeCoordFromMappedCoord(quadraturePoints[iQuad]);

      bool isAccepted = (acceptedData[idx] >= 0.);
      if(isAccepted)
      {
        const RealVector& otherState = originalData[acceptedIdx];
        const RealVector& pastTOtherState = pastData[acceptedIdx];
        tOtherState = varTransfo->transform(faces, coord, otherState,pastTOtherState);
        interfaceData[acceptedIdx] = *tOtherState;
        acceptedIdx++;
      }
      idx++;
    }
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::transformReceivedStatesData()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  ///Get Received DataHandle
  const std::string socketDataName = getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName,"States");
  DataHandle< RealVector> interfaceData =
      _sockets.getSocketSink<RealVector>(socketDataName)->getDataHandle();
  DataHandle< RealVector> originalData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_ORIGINAL")->getDataHandle();
  DataHandle< RealVector> pastData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_PAST")->getDataHandle();

  const std::string socketAcceptedName = getMethodData().getThisCoupledAcceptedName(_interfaceName,_currentTrsName,"States");
  DataHandle< CFreal> acceptedData =
      _sockets.getSocketSink< CFreal>(socketAcceptedName)->getDataHandle();
  const std::string socketConnectName = getMethodData().getThisCoupledConnectedName(_interfaceName,_currentTrsName,"States");
  DataHandle< vector<SubSysCouplerData::GeoEntityIdx> > connectData =
      _sockets.getSocketSink< vector<SubSysCouplerData::GeoEntityIdx> >(socketConnectName)->getDataHandle();

  Common::SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_currentTrsName);
  Common::SafePtr< vector<CFuint> > const trsStates = trs->getStatesInTrs();
  std::vector<State*>::iterator itd;

  /// Transform before setting the value to the DataHandle
  /// Call a variableTransformer
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  CFuint idx = 0;
  CFuint acceptedIdx = 0;

  RealVector* tOtherState(CFNULL);
  RealVector coord(PhysicalModelStack::getActive()->getDim());

  for (CFuint iState = 0; iState < trsStates->size(); ++iState)
  {
    State *const currState = states[(*trsStates)[iState]];
    coord = currState->getCoordinates();

    bool isAccepted = (acceptedData[idx] >= 0.);
    if(isAccepted)
    {
      const RealVector& otherState = originalData[acceptedIdx];
      const RealVector& pastTOtherState = pastData[acceptedIdx];
      tOtherState = varTransfo->transform(connectData[idx], coord, *currState, otherState, pastTOtherState);
      interfaceData[acceptedIdx] = (*tOtherState);

      acceptedIdx++;
    }
    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::transformReceivedNodalData()
{
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  ///Get Received DataHandle
  const std::string socketDataName = getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName,"Nodes");
  DataHandle< RealVector> interfaceData =
      _sockets.getSocketSink<RealVector>(socketDataName)->getDataHandle();
  DataHandle< RealVector> originalData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_ORIGINAL")->getDataHandle();
  DataHandle< RealVector> pastData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_PAST")->getDataHandle();

  const std::string socketAcceptedName = getMethodData().getThisCoupledAcceptedName(_interfaceName,_currentTrsName,"Nodes");
  DataHandle< CFreal> acceptedData =
      _sockets.getSocketSink< CFreal>(socketAcceptedName)->getDataHandle();
  const std::string socketConnectName = getMethodData().getThisCoupledConnectedName(_interfaceName,_currentTrsName,"Nodes");
  DataHandle< vector<SubSysCouplerData::GeoEntityIdx> > connectData =
      _sockets.getSocketSink< vector<SubSysCouplerData::GeoEntityIdx> >(socketConnectName)->getDataHandle();

  Common::SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_currentTrsName);
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();
  std::vector<Node*>::iterator itd;

  /// Transform before setting the value to the DataHandle
  /// Call a variableTransformer
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  CFuint idx = 0;
  CFuint acceptedIdx = 0;

  RealVector* tOtherState(CFNULL);
  RealVector coord(PhysicalModelStack::getActive()->getDim());

  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node *const currNode = nodes[(*trsNodes)[iNode]];

    bool isAccepted = (acceptedData[idx] >= 0.);
    if(isAccepted)
    {
      const RealVector& otherState = originalData[acceptedIdx];
      const RealVector& pastTOtherState = pastData[acceptedIdx];
      tOtherState = varTransfo->transform(connectData[idx], *currNode, otherState, pastTOtherState);
      interfaceData[acceptedIdx] = (*tOtherState);

      acceptedIdx++;
    }
    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::readFromDataHandle(const std::string dataHandleName)
{
  CFAUTOTRACE;

  DataHandle< RealVector> interfaceData =
    _sockets.getSocketSink<RealVector>(dataHandleName)->getDataHandle();

  DataHandle< RealVector> interfacePastData =
    _sockets.getSocketSink<RealVector>(dataHandleName + "_PAST")->getDataHandle();

  DataHandle< RealVector> originalData =
    _sockets.getSocketSink<RealVector>(dataHandleName + "_ORIGINAL")->getDataHandle();

  const std::string otherNamespace = getMethodData().getCoupledNameSpaceName(_interfaceName);
  const std::string datahandleName = otherNamespace + "_" + dataHandleName;

  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance().getNamespace(otherNamespace);
  Common::SafePtr<MeshData> otherMeshData = MeshDataStack::getInstance().getEntryByNamespace(otherNsp);

  DataHandle< RealVector> otherInterfaceData =
       otherMeshData->getDataStorage()->getData<RealVector>(datahandleName);

  cf_assert(otherInterfaceData.size() == interfaceData.size());
  cf_assert(otherInterfaceData.size() == interfacePastData.size());

  const CFuint nbStates = interfaceData.size();
  for (CFuint iState=0; iState < nbStates; ++iState)
  {
    //First backup past Data
    if(SubSystemStatusStack::getActive()->isSubIterationFirstStep())
    {
      interfacePastData[iState] = interfaceData[iState];
    }

    //Read the new data
    originalData[iState] = otherInterfaceData[iState];
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::readFile(const std::string dataFileName, const std::string acceptedFileName, const std::string dataHandleName, const std::string acceptedDataHandleName)
{
  CFAUTOTRACE;

  DataHandle< CFuint> parallelDataIndex =
    _sockets.getSocketSink<CFuint>(acceptedDataHandleName + "PAR")->getDataHandle();

  DataHandle< RealVector> interfaceData =
    _sockets.getSocketSink<RealVector>(dataHandleName)->getDataHandle();

  DataHandle< RealVector> interfacePastData =
    _sockets.getSocketSink<RealVector>(dataHandleName + "_PAST")->getDataHandle();

  DataHandle< RealVector> originalData =
    _sockets.getSocketSink<RealVector>(dataHandleName + "_ORIGINAL")->getDataHandle();

  //Read Accepted File
  using namespace boost::filesystem;
  path fname = Environment::DirPaths::getInstance().getResultsDir() / path ( acceptedFileName );

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(fname);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;
  std::vector<CFreal> tempIsAccepted(0);

  // Read nb states
  getWordsFromLine(fin,line,lineNb,words);
  cf_assert(words.size() == 1);
  CFuint nbStates = Common::StringOps::from_str<CFint>(words[0]);
  tempIsAccepted.resize(nbStates);
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    getWordsFromLine(fin,line,lineNb,words);
    tempIsAccepted[iState] = Common::StringOps::from_str<CFreal>(words[0]);
  }
  fhandle->close();



  //Read DataHandle
  path dataFname = Environment::DirPaths::getInstance().getResultsDir() / path ( dataFileName );

  Common::SelfRegistPtr<Environment::FileHandlerInput> dataFhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& dataFin = fhandle->open(dataFname);

  lineNb = 0;

  // Read nb states
  getWordsFromLine(dataFin,line,lineNb,words);
  cf_assert(words.size() == 1);
  nbStates = Common::StringOps::from_str<CFint>(words[0]);

  //Read the new data
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {

    if(tempIsAccepted[iState] >= 0.){
      getWordsFromLine(dataFin,line,lineNb,words);

      //Only store the transfered value if the data
      //comes from the processor who accepted the data
      if(parallelDataIndex[iState] == _iProc){
        //First backup past Data
        if(SubSystemStatusStack::getActive()->isSubIterationFirstStep())
        {
          interfacePastData[iState] = interfaceData[iState];
        }

        cf_assert(words.size() == (originalData[iState]).size());
        for (CFuint j=0; j<words.size();++j)
        {
          (originalData[iState])[j] =  Common::StringOps::from_str<CFreal>(words[j]);
        }
      }
    }
  }

  dataFhandle->close();
}


//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::prepareNormFile(const std::string dataHandleName)
{
  CFAUTOTRACE;

  DataHandle< RealVector> interfaceData =
    _sockets.getSocketSink<RealVector>(dataHandleName)->getDataHandle();

  const CFuint transferedSize = getMethodData().getThisTransferedSize(_interfaceName);
  const CFuint nbVars =
    getMethodData().getPostVariableTransformer(_interfaceName)->getTransformedSize(transferedSize);

  SimulationStatus::getInstance().addCouplingResidualNames(dataHandleName);
  CFout <<"Adding residual name: " << dataHandleName <<"\n";

  if (PE::GetPE().GetRank() == 0) {

    boost::filesystem::path fname =
      Environment::DirPaths::getInstance().getResultsDir() / ("Convergence_" + dataHandleName + ".plt");

    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(fname);

    fout << "TITLE = Coupling Convergence History" << "\n";
    fout << "VARIABLES = ";
    fout << "Iteration ";
    for(CFuint iVar = 0; iVar < nbVars ; iVar++)
    {
      fout << "L2Norm_Var[" <<iVar<<"] ";
    }
    for(CFuint iVar = 0; iVar < nbVars ; iVar++)
    {
      fout << "LInfNorm_Var[" <<iVar<<"] ";
    }
    fout << "\n";
    fhandle->close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::computeNorm(const std::string dataHandleName)
{
  CFAUTOTRACE;

  DataHandle< RealVector> interfaceData =
    _sockets.getSocketSink<RealVector>(dataHandleName)->getDataHandle();

  DataHandle< RealVector> interfacePastData =
    _sockets.getSocketSink<RealVector>(dataHandleName + "_PAST")->getDataHandle();

  cf_assert(interfaceData.size() == interfacePastData.size());

  // Read nb states
  const CFuint nbStates = interfaceData.size();

  cf_assert(nbStates > 0);
  const CFuint nbVars = interfaceData[0].size();
  RealVector delta(nbVars);
  delta = 0.;

  //Compute the square of the difference (for the L2 Norm)
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    cf_assert(interfaceData[iState].size() == nbVars);
    cf_assert(interfacePastData[iState].size() == nbVars);
    for (CFuint iVar = 0; iVar < nbVars; ++iVar)
    {
      const CFreal diff = std::fabs(interfaceData[iState][iVar] - interfacePastData[iState][iVar]);
      delta[iVar] += diff*diff;
    }
  }

  //Compute the maximum value of delta (for the LInf Norm)
  RealVector maxDelta(nbVars);
  maxDelta = 0.;
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    for (CFuint iVar = 0; iVar < nbVars; ++iVar)
    {
      const CFreal diff = std::fabs(interfaceData[iState][iVar] - interfacePastData[iState][iVar]);
      maxDelta[iVar] = std::max(diff,maxDelta[iVar]);
    }
  }

  /// @TODO temporary solution ....
  RealVector totalMaxDelta(nbVars);
  RealVector totalDelta(nbVars);
  for (CFuint iVar = 0; iVar < nbVars; ++iVar)
  {
    CFdouble tempDelta = 0.0;
    CFdouble tempMaxDelta = 0.0;
    if (PE::GetPE().GetProcessorCount() > 1) {
      CFdouble actualMaxDelta  = maxDelta[iVar];
      CFdouble actualDelta     = delta[iVar];
#ifdef CF_HAVE_MPI
      MPI_Allreduce(&actualDelta, &tempDelta, 1, MPI_DOUBLE, MPI_SUM,
		  PE::GetPE().GetCommunicator());

      MPI_Allreduce(&actualMaxDelta, &tempMaxDelta, 1, MPI_DOUBLE, MPI_MAX,
		  PE::GetPE().GetCommunicator());
#endif
    }
    else {
      tempMaxDelta  = maxDelta[iVar];
      tempDelta     = delta[iVar];
    }
    totalMaxDelta[iVar] = tempMaxDelta;
    totalDelta[iVar] = tempDelta;
  }


  CFreal L2norm = 0.;
  if (PE::GetPE().GetRank() == 0) {

    boost::filesystem::path fname =
      Environment::DirPaths::getInstance().getResultsDir() / ("Convergence_" + dataHandleName + ".plt");

    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(fname, ios::app);

    const CFuint iteration = SubSystemStatusStack::getActive()->getNbIter();

    fout << iteration;
    for (CFuint iVar = 0; iVar < nbVars; ++iVar){
      L2norm = log10(sqrt(totalDelta[iVar]));
      if(totalDelta[iVar] == 0.) L2norm = 0.;
  //    CFout << "Interface: " << _interfaceName << " - L2 Norm for variable["<< iVar << "]: " << L2norm << "\n";
      fout << " " << L2norm;
    }
    for (CFuint iVar = 0; iVar < nbVars; ++iVar){
      const CFreal LInfNorm = totalMaxDelta[iVar];
  //    CFout << "Interface: " << _interfaceName << " - LInf Norm for variable["<< iVar << "]: " << LInfNorm << "\n";
      fout << " " <<LInfNorm;
    }
    fout << "\n";
    fhandle->close();
  }

  ///@todo here should pass a vector!!
  SimulationStatus::getInstance().setCouplingResidual(L2norm,dataHandleName);
}

//////////////////////////////////////////////////////////////////////////////

void StdReadDataTransfer::getWordsFromLine(ifstream& fin,
                                      std::string& line,
                                      CFuint&  lineNb,
                                      vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

  } // namespace Numerics

} // namespace COOLFluiD
