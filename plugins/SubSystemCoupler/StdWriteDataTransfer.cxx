
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
#include "SubSystemCoupler/StdWriteDataTransfer.hh"

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

MethodCommandProvider<StdWriteDataTransfer, SubSysCouplerData, SubSystemCouplerModule> StdWriteDataTransferProvider("StdWriteDataTransfer");

//////////////////////////////////////////////////////////////////////////////

StdWriteDataTransfer::StdWriteDataTransfer(const std::string& name) :
  CouplerCom(name),
  _sockets(),
  socket_states("states"),
  socket_pastStates("pastStates")
{
}

//////////////////////////////////////////////////////////////////////////////

StdWriteDataTransfer::~StdWriteDataTransfer()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdWriteDataTransfer::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();;

  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdWriteDataTransfer::configure ( Config::ConfigArgs& args )
{
  CouplerCom::configure(args);

  typedef std::vector<Common::SafePtr<CommandGroup> > InterfaceList;
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
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdWriteDataTransfer::execute()
{
  CFAUTOTRACE;
  
  const std::string nsp = getMethodData().getNamespace();
  for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++iProc) {
    executeWrite(iProc);
  }

}

//////////////////////////////////////////////////////////////////////////////

void StdWriteDataTransfer::executeWrite(const CFuint iProc)
{

  CFAUTOTRACE;

  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();

  const std::string interfaceName = getCommandGroupName();
  const vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);
  const vector<std::string> coordTypes = getMethodData().getOtherCoordType(interfaceName);
  Common::SafePtr<PreVariableTransformer> varTransfo = getMethodData().getPreVariableTransformer(interfaceName);

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealVector otherState(nbEqs);
  RealVector* tOtherState(CFNULL);
  vector<SubSysCouplerData::GeoEntityIdx> faces(1);

  for (CFuint iTRS=0; iTRS < otherTrsNames.size(); iTRS++)
  {
    //Computing the values to be transfered to the other subsystem
    const vector<std::string> socketDataNames = getMethodData().getOtherCoupledDataName(interfaceName, otherTrsNames[iTRS],iProc);
    cf_assert(coordTypes.size() == socketDataNames.size());
    for(CFuint iType=0; iType < socketDataNames.size();iType++)
    {
      CFLogDebugMin("Writing data: " << socketDataNames[iType] << "\n");

      SubSysCouplerData::CoupledGeoEntities* coupledGeoEntities =
        getMethodData().getCoupledInterfaces(interfaceName,iTRS, iType, iProc);

      /// Get Datahandle of states of the Interface
      DataHandle< RealVector> interfaceData =
        _sockets.getSocketSink<RealVector>(socketDataNames[iType])->getDataHandle();

      const CFuint otherNbStates = interfaceData.size();
      varTransfo->setCurrentInterface(interfaceName, otherTrsNames[iTRS], coordTypes[iType]);
      varTransfo->setNbStates(otherNbStates);

      /// Get the geometric entity builder
      Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
      geoBuilder = getMethodData().getStdTrsGeoBuilder();

      StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

      for (CFuint iState = 0; iState < otherNbStates; ++iState)
      {
        RealVector& shapeFunctions = (*coupledGeoEntities)[iState].third;
        RealVector& coord = (*coupledGeoEntities)[iState].fourth;

        ///Create the Geometric Entity to get the face corresponding to the idx/TRS given in .first and .second
        // build the GeometricEntity
        geoData.idx = (*coupledGeoEntities)[iState].second;
        geoData.trs = (*coupledGeoEntities)[iState].first;
        GeometricEntity& currFace = *geoBuilder->buildGE();

        faces[0].first = (*coupledGeoEntities)[iState].first;
        faces[0].second = geoData.idx;

        const CFuint nbNodes = shapeFunctions.size();

        std::vector<State*>* states = currFace.getStates();

        cf_assert(nbNodes == states->size());
        cf_assert ((*(*states)[0]).size() == nbEqs);
        for (CFuint j = 0; j < nbEqs; ++j) {
          otherState[j] = shapeFunctions[0] * ((*(*states)[0])[j]);
          for (CFuint k = 1; k < nbNodes; ++k) {
            otherState[j] += shapeFunctions[k] * ((*((*states)[k]))[j]);
          }
        }

        ///release GeometricEntity
        geoBuilder->releaseGE();

        varTransfo->setStateIndex(iState);
        tOtherState = varTransfo->preTransform(faces, coord, otherState, shapeFunctions);
        interfaceData[iState] = (*tOtherState);

        CFLogDebugMin("=============================================================" << "\n");
        CFLogDebugMin("=========== " << socketDataNames[iType] <<" ====================================" << "\n");
        CFLogDebugMin("================================================================" << "\n");
        CFLogDebugMin("Writing node: " << coord << "\n");
        CFLogDebugMin("...iState: " << iState << "\n");
        CFLogDebugMin("...socketName: " << socketDataNames[iType] << "\n");
        CFLogDebugMin("...ShapeFunction: " << shapeFunctions << "\n");
        CFLogDebugMin("...nbStates: " << states->size() << "\n");
        CFLogDebugMin("...State[0] coord: " << ((*states)[0])->getCoordinates() << "\n");
        CFLogDebugMin("...State[1] coord: " << ((*states)[1])->getCoordinates() << "\n");
        CFLogDebugMin("...State[0]: " << *((*states)[0]) << "\n");
        CFLogDebugMin("...State[1]: " << *((*states)[1]) << "\n");
//        CFLogDebugMin("...State[2]: " << *((*states)[2]) << "\n");
        CFLogDebugMin("...Non Transformed State: " << otherState << "\n");
        CFLogDebugMin("...Transformed State: " << tOtherState << "\n");

      }

      ///Write the datahandle to a file
      if(getMethodData().isTransferFiles()) writeFile(socketDataNames[iType]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdWriteDataTransfer::writeFile(const std::string socketName)
{
  CFAUTOTRACE;

  DataHandle< RealVector> interfaceData =
    _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();

  boost::filesystem::path nameOutputFile =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(socketName);

// CFout << "Writing data file: " <<socketName <<"\n";
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(nameOutputFile);

  const CFuint nbStates = interfaceData.size();

  fout << nbStates << "\n";
  for (CFuint i = 0; i < nbStates; ++i) {
    fout << interfaceData[i] << "\n";
  }

  fhandle->close();
// CFout << "Done writing data file: " <<socketName <<"\n";
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

  } // namespace Numerics

} // namespace COOLFluiD
