
#include <boost/filesystem/path.hpp>

#include "Common/FilesystemException.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodCommand.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/VolumeIntegratorImpl.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/FVMCCReadDataTransfer.hh"

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

MethodCommandProvider<FVMCCReadDataTransfer, SubSysCouplerData, SubSystemCouplerModule> FVMCCReadDataTransferProvider("FVMCCReadDataTransfer");

//////////////////////////////////////////////////////////////////////////////

FVMCCReadDataTransfer::FVMCCReadDataTransfer(const std::string& name) :
  StdReadDataTransfer(name),
  socket_nstates("nstates"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCCReadDataTransfer::~FVMCCReadDataTransfer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FVMCCReadDataTransfer::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdReadDataTransfer::needsSockets();

  result.push_back(&socket_nstates);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCReadDataTransfer::transformReceivedGhostData()
{
  CFAUTOTRACE;

  //Get Received DataHandle
  const std::string socketDataName = getMethodData().getThisCoupledDataName(_interfaceName,_currentTrsName,"Ghost");
  DataHandle< RealVector> interfaceData =
      _sockets.getSocketSink<RealVector>(socketDataName)->getDataHandle();
  DataHandle< RealVector> originalData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_ORIGINAL")->getDataHandle();
  DataHandle< RealVector> pastData =
      _sockets.getSocketSink<RealVector>(socketDataName + "_PAST")->getDataHandle();

  const std::string socketAcceptedName = getMethodData().getThisCoupledAcceptedName(_interfaceName,_currentTrsName,"Ghost");
  DataHandle< CFreal> acceptedData =
      _sockets.getSocketSink< CFreal>(socketAcceptedName)->getDataHandle();

  // Transform before setting the value to the DataHandle
  // Call a variableTransformer
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  // Get the geometric entity builder
  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  Common::SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(_currentTrsName);
  geoData.trs = trs;
  geoData.isBFace = true;
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  RealVector* tOtherState(CFNULL);
  vector<SubSysCouplerData::GeoEntityIdx> faces(1);
  RealVector coord(PhysicalModelStack::getActive()->getDim());

  // Fill the datahandle
  // index going through the whole TRS
  CFuint idx = 0;
  // index going through the accepted points in the TRS
  CFuint acceptedIdx = 0;
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
  {
    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity* currFace = geoBuilder->buildGE();

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
        tOtherState = varTransfo->transform(faces, coord, otherState, pastTOtherState);
        interfaceData[acceptedIdx] = *tOtherState;

        acceptedIdx++;
      }
      idx++;
    }
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCReadDataTransfer::transformReceivedNodalData()
{
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
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
  std::vector<State*>::iterator itd;

  /// Transform before setting the value to the DataHandle
  /// Call a variableTransformer
  Common::SafePtr<PostVariableTransformer> varTransfo = getMethodData().getPostVariableTransformer(_interfaceName);

  RealVector* tOtherState(CFNULL);
  RealVector coord(PhysicalModelStack::getActive()->getDim());

  CFuint idx = 0;
  CFuint acceptedIdx = 0;
  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node *const currNode = nodes[(*trsNodes)[iNode]];
    const RealVector& currState = nstates[(*trsNodes)[iNode]];
    coord = (*currNode);

    bool isAccepted = (acceptedData[idx] >= 0.);
    if(isAccepted)
    {
      const RealVector& otherState = originalData[acceptedIdx];
      const RealVector& pastTOtherState = pastData[acceptedIdx];
      tOtherState = varTransfo->transform(connectData[idx], coord, currState, otherState, pastTOtherState);
      interfaceData[acceptedIdx] = *tOtherState;

      acceptedIdx++;
    }
    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

  } // namespace Numerics

} // namespace COOLFluiD
