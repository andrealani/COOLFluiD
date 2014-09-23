
#include "Framework/MethodCommandProvider.hh"
#include "PLaS/PLaSModule.hh"
#include "PLaS/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSetup,PLaSTrackingData,PLaSModule > stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFLog(INFO,"PLaSTracking: initialize library...\n");
  getMethodData().PLaS_Init();
  CFLog(VERBOSE,"PLaSTracking: initialize library.\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

