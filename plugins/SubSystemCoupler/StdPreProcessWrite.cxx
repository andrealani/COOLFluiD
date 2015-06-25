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
#include "SubSystemCoupler/StdPreProcessWrite.hh"
#include "Common/BadValueException.hh"

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

MethodCommandProvider<StdPreProcessWrite, SubSysCouplerData, SubSystemCouplerModule> stdPreProcessWriteProvider("StdPreProcessWrite");

//////////////////////////////////////////////////////////////////////////////

StdPreProcessWrite::StdPreProcessWrite(const std::string& name) : CouplerCom(name),
  _sockets(),
  _nbCoupledStates(0),
  _alreadyReadSockets(0),
  socket_states("states"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdPreProcessWrite::~StdPreProcessWrite()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

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
      for(; trs != trsNames.end(); ++trs)
      {
        const std::string trsName = *trs;

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
          _sockets.createSocketSource<CFreal>(socketAcceptNames[iType]);
          _sockets.createSocketSource<RealVector>(socketDataNames[iType]);
          _sockets.createSocketSource<RealVector>(socketDataNames[iType] + "_PAST");
          _sockets.createSocketSource<RealVector>(socketDataNames[iType] + "_ORIGINAL");
          _sockets.createSocketSource<RealVector>(socketCoordNames[iType]);
          _sockets.createSocketSource<std::vector<SubSysCouplerData::GeoEntityIdx> >(socketConnectNames[iType]);
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

void StdPreProcessWrite::executeOnTrs()
{
  CFAUTOTRACE;
  // CFout << "Execute PreProcessWrite BEGIN\n";
  
  const std::string nsp = getMethodData().getNamespace();
  for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++iProc) {
    executeWriteOnTrs(iProc);
  }
// CFout << "Execute PreProcessWrite END\n";
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::executeWriteOnTrs(const CFuint iProc)
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLog(VERBOSE,"StdPreProcessWrite::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // Names of the interfaces, subsystems
  const std::string interfaceName = getCommandGroupName();
  const std::string trsName = trs->getName();

  // socket for accepted data
  vector<std::string> socketAcceptNames = getMethodData().getThisCoupledAcceptedName(interfaceName,trsName);
  vector<std::string> socketDataNames = getMethodData().getThisCoupledDataName(interfaceName,trsName);
  vector<std::string> socketCoordNames = getMethodData().getThisCoupledCoordName(interfaceName,trsName);
  vector<std::string> socketConnectNames = getMethodData().getThisCoupledConnectedName(interfaceName,trsName);
  vector<std::string> coordTypes = getMethodData().getThisCoordType(interfaceName);

  for(CFuint iType=0;iType < socketAcceptNames.size();iType++)
  {
    DataHandle< RealVector> interfaceCoords =
      _sockets.getSocketSource<RealVector>(socketCoordNames[iType])->getDataHandle();

      if(coordTypes[iType] == "Nodes")
      {
        // resize nodal sockets
        _nbCoupledStates = getCurrentTRS()->getNodesInTrs()->size();
        interfaceCoords.resize(_nbCoupledStates);

        fillNodalCoordDataHandle(socketCoordNames[iType]);
        setNodeToFaceConnectivity(socketConnectNames[iType]);
      }
      if(coordTypes[iType] == "States")
      {
        // resize nodal sockets
        _nbCoupledStates = getCurrentTRS()->getStatesInTrs()->size();
        interfaceCoords.resize(_nbCoupledStates);

        fillStatesCoordDataHandle(socketCoordNames[iType]);
        setStateToFaceConnectivity(socketConnectNames[iType]);
      }

      if(coordTypes[iType] == "Gauss")
      {
        setQuadratureDataHandles();

        interfaceCoords.resize(_nbCoupledStates);
        fillQuadratureCoordDataHandle(socketCoordNames[iType]);
      }
      if(coordTypes[iType] == "Ghost")
      {
        _nbCoupledStates = getCurrentTRS()->getLocalNbGeoEnts();
        interfaceCoords.resize(_nbCoupledStates);

        fillGhostCoordDataHandle(socketCoordNames[iType]);
      }

      if(getMethodData().isInterfaceGeometryNonMatching(interfaceName))
      {
        CFuint dim = 0;
        if(interfaceCoords.size() > 0) dim = interfaceCoords[0].size();
        const CFreal rotationAngle = getMethodData().getInterfaceRotation(interfaceName);
        const RealVector& translationVector = getMethodData().getInterfaceTranslation(interfaceName);
        RealVector coord(dim);
    //  const vector<CFreal> rotationCenter(dim);
    //  rotationCenter = 0.;
        CFreal cosT = cos(rotationAngle);
        CFreal sinT = sin(rotationAngle);

        if(rotationAngle != 0) {cf_assert(dim ==2);}

        for (CFuint idx = 0; idx< interfaceCoords.size(); ++idx) {
          // Modify the coordinates: Rotation (2D)
          ///@todo extend the rotation for 3D
          coord[0] = cosT*(interfaceCoords[idx])[0] + sinT*(interfaceCoords[idx])[1];
          coord[1] = -sinT*(interfaceCoords[idx])[0] + cosT*(interfaceCoords[idx])[1];

          for(CFuint iCoord=0;iCoord<interfaceCoords[idx].size();++iCoord)
          {
            (interfaceCoords[idx])[iCoord] = coord[iCoord] + translationVector[iCoord];
          }
        }//end for
      }//end if
      writeFile(socketCoordNames[iType]);

    // Resize the datahandles, now that we know the number of coupled states
/*    DataHandle<CFreal> isAccepted =
      _sockets.getSocketSource<CFreal>(socketAcceptNames[iType])->getDataHandle();
    DataHandle< RealVector> interfaceData =
      _sockets.getSocketSource<RealVector>(socketDataNames[iType])->getDataHandle();
    DataHandle< RealVector> interfacePastData =
      _sockets.getSocketSource<RealVector>(socketDataNames[iType] + "_PAST")->getDataHandle();
    DataHandle< RealVector> interfaceOrigData =
      _sockets.getSocketSource<RealVector>(socketDataNames[iType] + "_ORIGINAL")->getDataHandle();

    isAccepted.resize(_nbCoupledStates);
    interfaceData.resize(_nbCoupledStates);
    interfacePastData.resize(_nbCoupledStates);
    interfaceOrigData.resize(_nbCoupledStates);

    //Resize the data to the correct size
    cf_assert(interfaceData.size() == interfacePastData.size());
    cf_assert(interfaceData.size() == interfaceOrigData.size());

    const CFuint transferedSize = getMethodData().getThisTransferedSize(interfaceName);
    const CFuint transformedSize =
      getMethodData().getPostVariableTransformer(interfaceName)->getTransformedSize(transferedSize);

    for(CFuint iData=0;iData< interfaceData.size(); iData++)
    {
      interfaceOrigData[iData].resize(transferedSize);
      interfaceData[iData].resize(transformedSize);
      interfacePastData[iData].resize(transformedSize);
    }

    //Reset the data and pastData to zero
    cf_assert(interfaceData.size() == interfacePastData.size());
    for(CFuint iData=0;iData< interfaceData.size(); iData++)
    {
      interfaceData[iData] = 0.;
      interfacePastData[iData] = 0.;
      interfaceOrigData[iData] = 0.;
    }*/
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::setQuadratureDataHandles()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  /// Compute size of datahandle
  CFuint totalSize = 0;
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {
    CFLogDebugMax("Face " << iGeoEnt << "\n");

    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity* currFace = geoBuilder->buildGE();

    // number of quadrature points associated with this face
    const CFuint nbPoints =
      getMethodData().getCollaborator<SpaceMethod>()->getVolumeIntegrator()->getSolutionIntegrator(currFace)->getIntegratorPattern().totalNbPts();

    totalSize += nbPoints;

    //release the GeometricEntity
    geoBuilder->releaseGE();
  }

  /// resize the data to be transfered TO THE CURRENT SubSystem
  _nbCoupledStates = totalSize;
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::fillStatesCoordDataHandle(const std::string& socketName)
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  /// Fill in the coordinates datahandle
  DataHandle<RealVector> interfaceCoord =
    _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsStates = trs->getStatesInTrs();

  CFuint idx = 0;
  for (CFuint iState = 0; iState < trsStates->size(); ++iState)
  {
    State *const currState = states[(*trsStates)[iState]];

    interfaceCoord[idx].resize(currState->getCoordinates().size());
    interfaceCoord[idx] = currState->getCoordinates();
    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::setStateToFaceConnectivity(const std::string& socketName)
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsStates = trs->getStatesInTrs();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  /// resize the socket of pointers to GeomEntity
  DataHandle<std::vector<SubSysCouplerData::GeoEntityIdx> > interfaceStateToFaceConn =
    _sockets.getSocketSource<std::vector<SubSysCouplerData::GeoEntityIdx> >(socketName)->getDataHandle();
  interfaceStateToFaceConn.resize(trsStates->size());

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  geoData.trs = trs;
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  CFuint idx = 0;
  for (CFuint iState = 0; iState < trsStates->size(); ++iState)
  {
    State *const currState = states[(*trsStates)[iState]];
    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {
      // build the GeometricEntity
      geoData.idx = iGeoEnt;

      GeometricEntity* currFace = geoBuilder->buildGE();
      bool contained = currFace->containState(currState);
      if(contained)
      {
        SubSysCouplerData::GeoEntityIdx geoEntIdx;
        geoEntIdx.first = trs;
        geoEntIdx.second = iGeoEnt;
        interfaceStateToFaceConn[idx].push_back(geoEntIdx);
      }
    //release the GeometricEntity
    geoBuilder->releaseGE();
    }

    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::fillNodalCoordDataHandle(const std::string& socketName)
{
  CFAUTOTRACE;

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  /// Fill in the coordinates datahandle
  DataHandle<RealVector> interfaceCoord =
    _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();

  CFuint idx = 0;
  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node *const currNode = nodes[(*trsNodes)[iNode]];

    interfaceCoord[idx].resize(currNode->size());
    interfaceCoord[idx] = *currNode;
    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::setNodeToFaceConnectivity(const std::string& socketName)
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  /// resize the socket of pointers to GeomEntity
  DataHandle<std::vector<SubSysCouplerData::GeoEntityIdx> > interfaceNodeToFaceConn =
    _sockets.getSocketSource<std::vector<SubSysCouplerData::GeoEntityIdx> >(socketName)->getDataHandle();
  interfaceNodeToFaceConn.resize(trsNodes->size());

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  geoData.trs = trs;
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  CFuint idx = 0;
  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node *const currNode = nodes[(*trsNodes)[iNode]];
    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {
      // build the GeometricEntity
      geoData.idx = iGeoEnt;

      GeometricEntity* currFace = geoBuilder->buildGE();
      bool contained = currFace->containNode(currNode);
      if(contained)
      {
        SubSysCouplerData::GeoEntityIdx geoEntIdx;
        geoEntIdx.first = trs;
        geoEntIdx.second = iGeoEnt;
        interfaceNodeToFaceConn[idx].push_back(geoEntIdx);
      }
    //release the GeometricEntity
    geoBuilder->releaseGE();
    }

    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::fillQuadratureCoordDataHandle(const std::string& socketName)
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  /// Get the coordinates datahandle
  DataHandle<RealVector> interfaceCoord =
    _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  geoData.trs = trs;
  const CFuint nbGeos = trs->getLocalNbGeoEnts();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  RealVector quadCoord(nbDim);

  /// Fill the datahandle
  CFuint idx = 0;

    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {
      // build the GeometricEntity
      geoData.idx = iGeoEnt;
      GeometricEntity* currFace = geoBuilder->buildGE();

      // gets the correct integrator corresponding to the
      // shape of the element where you want to integrate
      VolumeIntegratorImpl *const impl =
        getMethodData().getCollaborator<SpaceMethod>()->getVolumeIntegrator()->getSolutionIntegrator(currFace);

      std::vector<RealVector> quadraturePoints = impl->getQuadraturePointsCoordinates();

      for (CFuint iQuad = 0; iQuad < quadraturePoints.size(); ++iQuad)
      {
//      CFout << "QuadPoint: " << quadraturePoints[iQuad] << "\n";
        quadCoord = currFace->computeCoordFromMappedCoord(quadraturePoints[iQuad]);
//        CFout <<  "QuadPoint Coord: " << quadCoord << "\n";
        interfaceCoord[idx].resize(quadCoord.size());
        interfaceCoord[idx] = quadCoord;
        idx++ ;
      }
    //release the GeometricEntity
    geoBuilder->releaseGE();
    }
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::writeFile(const std::string& socketName)
{
  CFAUTOTRACE;

  DataHandle<RealVector> interfaceCoord =
    _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

  const std::string nsp = getMethodData().getNamespace();
  std::string fileName = socketName;
  const bool isParallel = Common::PE::GetPE().IsParallel();
  if(isParallel) fileName += ".P" + StringOps::to_str(Common::PE::GetPE().GetRank(nsp));

  boost::filesystem::path nameOutputFile =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(fileName);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(nameOutputFile);

  const CFuint nbStates = interfaceCoord.size();

  fout << nbStates << "\n";
  for (CFuint i = 0; i < nbStates; ++i) {
    fout << interfaceCoord[i] << "\n";
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void StdPreProcessWrite::getWordsFromLine(ifstream& fin,
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
StdPreProcessWrite::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = _sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdPreProcessWrite::needsSockets()
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
