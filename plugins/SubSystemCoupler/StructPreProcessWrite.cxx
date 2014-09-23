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
#include "SubSystemCoupler/StructPreProcessWrite.hh"
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

MethodCommandProvider<StructPreProcessWrite, SubSysCouplerData, SubSystemCouplerModule> StructPreProcessWriteProvider("StructPreProcessWrite");

//////////////////////////////////////////////////////////////////////////////

StructPreProcessWrite::StructPreProcessWrite(const std::string& name) :
  StdPreProcessWrite(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StructPreProcessWrite::~StructPreProcessWrite()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructPreProcessWrite::fillNodalCoordDataHandle(const std::string& socketName)
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

    cf_assert(currState->size() == currState->getCoordinates().size());
    interfaceCoord[idx] = currState->getCoordinates() + (*currState);

    idx++ ;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD
