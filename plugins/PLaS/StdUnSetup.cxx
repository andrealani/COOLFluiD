
#include "Framework/MethodCommandProvider.hh"
#include "PLaS/PLaSModule.hh"
#include "PLaS/StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdUnSetup,PLaSTrackingData,PLaSModule > stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFLog(INFO,"PLaSTracking: terminate library...\n");
  getMethodData().PLaS_Terminate();
  CFLog(VERBOSE,"PLaSTracking: terminate library.\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

