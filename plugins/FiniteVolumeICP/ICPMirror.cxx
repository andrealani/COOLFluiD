#include "FiniteVolumeICP/ICPMirror.hh"
#include "FiniteVolumeICP/FiniteVolumeICP.hh"
#include "FiniteVolume/MirrorVelocity.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ICPMirror<MirrorVelocity>, 
		      CellCenterFVMData, 
		      FiniteVolumeICPModule> 
mirrorICPFVMCCProvider("MirrorICPFVMCC");


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
