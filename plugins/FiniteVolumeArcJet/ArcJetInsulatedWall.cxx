#include "FiniteVolumeArcJet/ArcJetInsulatedWall.hh"
#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "ArcJet/ArcJetInductionTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::ArcJet;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ArcJetInsulatedWall<NoSlipWallIsothermalNSvt<ArcJetInductionTerm<EulerTerm> > >, 
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
noSlipWallIsothermalArcJetPvtFVMCCProvider("ArcJetInsulatedNoSlipWallFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
