#include "FiniteVolumeArcJet/ArcJetPhiInsulatedWall.hh"
#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVT.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVTLTE.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ArcJetPhiInsulatedWall<NoSlipWallIsothermalNSvt<MultiScalarTerm<EulerTerm> > >, 
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiInsulatedWallFVMCCProvider("ArcJetPhiInsulatedWallFVMCC");

MethodCommandProvider<ArcJetPhiInsulatedWall< SubInletEulerPvtVT<EulerVarSet> >,
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiInletFVMCC("ArcJetPhiInletFVMCC");

MethodCommandProvider<ArcJetPhiInsulatedWall< SubInletEulerPvtVTLTE >,
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiInletMassFlowFVMCC("ArcJetPhiInletMassFlowFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
