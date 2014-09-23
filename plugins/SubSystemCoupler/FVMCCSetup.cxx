#include "Framework/MethodCommandProvider.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "SubSystemCoupler/FVMCCSetup.hh"
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

MethodCommandProvider<FVMCCSetup, SubSysCouplerData, SubSystemCouplerModule> FVMCCSetupProvider("FVMCCSetup");

//////////////////////////////////////////////////////////////////////////////

FVMCCSetup::FVMCCSetup(const std::string& name) :
  StdSetup(name),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_gstates("gstates")
{
}

//////////////////////////////////////////////////////////////////////////////

FVMCCSetup::~FVMCCSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCCSetup::execute()
{
  CFAUTOTRACE;

  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();

  Common::SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FVMCCSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = StdSetup::needsSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD
