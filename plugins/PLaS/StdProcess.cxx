
#include "Framework/MethodCommandProvider.hh"
#include "PLaS/PLaSModule.hh"
#include "PLaS/StdProcess.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdProcess,PLaSTrackingData,PLaSModule > stdProcessProvider("StdProcess");

//////////////////////////////////////////////////////////////////////////////

void StdProcess::execute()
{
  CFLog(INFO,"PLaSTracking: run library...\n");
  getMethodData().PLaS_Run();
  CFLog(VERBOSE,"PLaSTracking: run library.\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

