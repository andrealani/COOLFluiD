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
#include "SubSystemCoupler/FVMCCPreProcessWrite.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
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

MethodCommandProvider<FVMCCPreProcessWrite, SubSysCouplerData, SubSystemCouplerModule> FVMCCPreProcessWriteProvider("FVMCCPreProcessWrite");

//////////////////////////////////////////////////////////////////////////////

FVMCCPreProcessWrite::FVMCCPreProcessWrite(const std::string& name) :
  StdPreProcessWrite(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCCPreProcessWrite::~FVMCCPreProcessWrite()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCPreProcessWrite::executeOnTrs()
{
  CFAUTOTRACE;

  StdPreProcessWrite::executeOnTrs();

}

//////////////////////////////////////////////////////////////////////////////

void FVMCCPreProcessWrite::fillGhostCoordDataHandle(const std::string& socketName)
{
  CFAUTOTRACE;

  /// Fill in the coordinates datahandle
  DataHandle<RealVector> interfaceCoord =
    _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

  RealVector bCoord(PhysicalModelStack::getActive()->getDim());

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  geoData.trs = getCurrentTRS();
  geoData.isBFace = true;

  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
  // unused //  const CFuint nbTrsStates = getCurrentTRS()->getStatesInTrs()->size();

  cf_assert(nbTrsFaces == interfaceCoord.size());

  //Loop over faces
  CFuint idx = 0;
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace)
  {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();

    State *const innerState = currFace->getState(0);
    State *const ghostState = currFace->getState(1);

    // coordinate of the boundary point
    bCoord = (innerState->getCoordinates() +
              ghostState->getCoordinates());
    bCoord *= 0.5;

    //Write the coordinates
    interfaceCoord[idx].resize(bCoord.size());
    interfaceCoord[idx] = bCoord;
    idx++ ;

    //release the face
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCPreProcessWrite::setNodeToFaceConnectivity(const std::string& socketName)
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  /// resize the socket of pointers to GeomEntity
  DataHandle<std::vector<SubSysCouplerData::GeoEntityIdx> > interfaceNodeToFaceConn =
    _sockets.getSocketSource<std::vector<SubSysCouplerData::GeoEntityIdx> >(socketName)->getDataHandle();
  interfaceNodeToFaceConn.resize(trsNodes->size());

  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  geoData.trs = getCurrentTRS();
  geoData.isBFace = true;

  const CFuint nbGeos = trs->getLocalNbGeoEnts();
  CFuint idx = 0;
  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node *const currNode = nodes[(*trsNodes)[iNode]];

    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {
      // build the GeometricEntity
      geoData.idx = iGeoEnt;

      GeometricEntity& currFace = *geoBuilder->buildGE();
      bool contained = currFace.containNode(currNode);
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

void FVMCCPreProcessWrite::fillNodalCoordDataHandle(const std::string& socketName)
{
  CFAUTOTRACE;

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // Fill in the coordinates datahandle
  DataHandle<RealVector> interfaceCoord =
    _sockets.getSocketSource<RealVector>(socketName)->getDataHandle();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = trs->getNodesInTrs();

  cf_assert(interfaceCoord.size() == trsNodes->size());
  CFuint idx = 0;
  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node *const currNode = nodes[(*trsNodes)[iNode]];
    interfaceCoord[idx].resize(currNode->size());
    interfaceCoord[idx] = (*currNode);
    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCPreProcessWrite::fillStatesCoordDataHandle(const std::string& socketName)
{
  CFAUTOTRACE;

  throw BadValueException (FromHere(),"Transfer at the states cannot be chosen for FVMCC method, you may want to switch to transfer at the nodes.");

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD
