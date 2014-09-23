
#include <boost/filesystem/path.hpp>

#include "Common/FilesystemException.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/VolumeIntegratorImpl.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/FVMCCWriteDataTransfer.hh"

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

MethodCommandProvider<FVMCCWriteDataTransfer, SubSysCouplerData, SubSystemCouplerModule> FVMCCWriteDataTransferProvider("FVMCCWriteDataTransfer");

//////////////////////////////////////////////////////////////////////////////

FVMCCWriteDataTransfer::FVMCCWriteDataTransfer(const std::string& name) :
  StdWriteDataTransfer(name),
  socket_nstates("nstates")
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCCWriteDataTransfer::~FVMCCWriteDataTransfer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FVMCCWriteDataTransfer::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdWriteDataTransfer::needsSockets();

  result.push_back(&socket_nstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCWriteDataTransfer::executeWrite(const CFuint iProc)
{

  CFAUTOTRACE;

  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();

  const std::string interfaceName = getCommandGroupName();
  ///@todo change with This TRSNames
  const vector<std::string> trsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);
  const vector<std::string> coordTypes = getMethodData().getThisCoordType(interfaceName);
  Common::SafePtr<PreVariableTransformer> varTransfo = getMethodData().getPreVariableTransformer(interfaceName);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealVector otherState(nbEqs);
  RealVector* tOtherState(CFNULL);
  vector<SubSysCouplerData::GeoEntityIdx> faces(1);
  for (CFuint iTRS=0; iTRS < trsNames.size(); iTRS++)
  {
    //Computing the values to be transfered to the other subsystem
    const vector<std::string> socketDataNames = getMethodData().getOtherCoupledDataName(interfaceName, trsNames[iTRS],iProc);
    for(CFuint iType=0; iType < socketDataNames.size();iType++)
    {
      SubSysCouplerData::CoupledGeoEntities* coupledGeoEntities =
        getMethodData().getCoupledInterfaces(interfaceName,iTRS, iType, iProc);

      /// Get Datahandle of states of the Interface
      DataHandle< RealVector> interfaceData =
        _sockets.getSocketSink<RealVector>(socketDataNames[iType])->getDataHandle();

      const CFuint otherNbStates = interfaceData.size();
      varTransfo->setCurrentInterface(interfaceName, trsNames[iTRS], coordTypes[iType]);
      varTransfo->setNbStates(otherNbStates);

      /// Get the geometric entity builder
      Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
      geoBuilder = getMethodData().getFaceTrsGeoBuilder();
      FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

      for (CFuint iState = 0; iState < otherNbStates; ++iState)
      {
        RealVector& shapeFunctions = (*coupledGeoEntities)[iState].third;
        RealVector& coord = (*coupledGeoEntities)[iState].fourth;
        // std::cout << "Writing node: " << coord << std::endl;
        // std::cout << "...i: " << i << std::endl;
        // std::cout << "...socketName: " << socketName << std::endl;
        // std::cout << "...ShapeFunction: " << shapeFunctions << std::endl;

        ///Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
        // build the GeometricEntity
        geoData.idx = (*coupledGeoEntities)[iState].second;
        geoData.trs = (*coupledGeoEntities)[iState].first;
        geoData.isBFace = true;
        GeometricEntity& currFace = *geoBuilder->buildGE();

        const std::vector<Node*>* nodes = currFace.getNodes();
        const CFuint nbNodes = nodes->size();
        std::vector<RealVector> nstates(nbNodes);

        //Get the nodal states cooresponding to the face nodes
        for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
          // std::cout << "iNode: " << iNode <<std::endl;
          // std::cout << "Node Local ID: " << (*nodes)[iNode]->getLocalID() <<std::endl;
          // std::cout << "Node: " << *((*nodes)[iNode]) <<std::endl;
          // std::cout << "Nodal State: " << nodalStates[(*nodes)[iNode]->getLocalID()] <<std::endl;
          nstates[iNode].resize(nbEqs);
          nstates[iNode] = nodalStates[(*nodes)[iNode]->getLocalID()];
        }

        //std::cout << "...nbStates: " << nstates.size() << std::endl;
        //std::cout << "...State[0]: " << nstates[0] << std::endl;
        //std::cout << "...State[1]: " << nstates[1] << std::endl;
        //std::cout << "...State[2]: " << nstates[2] << std::endl;
        faces[0].first = (*coupledGeoEntities)[iState].first;
        faces[0].second = geoData.idx;

        cf_assert(nbNodes == shapeFunctions.size());
        cf_assert ((nstates[0]).size() == nbEqs);

//      std::cout << "-------------------------------------" <<std::endl;
//      for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
//        std::cout << "Node: " << *((*nodes)[iNode]) <<std::endl;
//        std::cout << "NState: " << nstates[iNode] <<std::endl;
//        std::cout << "ShapeFunctions: " <<shapeFunctions[iNode] << std::endl;
//      }

        for (CFuint j = 0; j < nbEqs; ++j) {
          otherState[j] = shapeFunctions[0] * ((nstates[0])[j]);
          for (CFuint k = 1; k < nbNodes; ++k) {
            otherState[j] += shapeFunctions[k] * ((nstates[k])[j]);
          }
        }

        ///release GeometricEntity
        geoBuilder->releaseGE();

        //std::cout << "...Non Transformed State: " << otherState << std::endl;
        /// Transform before setting the value to the DataHandle
        varTransfo->setStateIndex(iState);
        tOtherState = varTransfo->preTransform(faces, coord, otherState, shapeFunctions);
        interfaceData[iState] = (*tOtherState);

        //std::cout << "...Transformed State: " << tOtherState << std::endl;
      }

      ///Write the datahandle to a file
      if(getMethodData().isTransferFiles()) writeFile(socketDataNames[iType]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

  } // namespace Numerics

} // namespace COOLFluiD
