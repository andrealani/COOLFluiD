#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolumeArcJet/ArcJetPhiOutlet.hh"
#include "FiniteVolumeArcJet/ArcJetPhiInsulatedWall.hh"
#include "FiniteVolumeArcJet/ArcJetPhiElectrode.hh"
#include "FiniteVolumeNavierStokes/SubBCTurb.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVTLTE.hh"
#include "FiniteVolumeTurb/NoSlipWallIsothermalNSKOmPvt.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"
#include "LTE/Euler3DPvtLTE.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

// outlet with given presssure
MethodCommandProvider<ArcJetPhiOutlet<SubBCTurb
		      <SubOutletEulerPvt<MultiScalarVarSet<Euler3DPvtLTE> >, 
		       MultiScalarVarSet<Euler3DPvtLTE>, false> >, 
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiOutlet3DTurbFVMCCProvider("ArcJetPhiOutlet3DTurbFVMCC");

// inlet with given mass flow
MethodCommandProvider<ArcJetPhiInsulatedWall<SubBCTurb
		      <SubInletEulerPvtVTLTE, MultiScalarVarSet<Euler3DPvtLTE>, true> >, 
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiInletTurbFVMCC("ArcJetPhiInletTurbFVMCC");

// isothermal wall
MethodCommandProvider<ArcJetPhiInsulatedWall<NoSlipWallIsothermalNSKOmPvt
		      <EulerVarSet, NavierStokesTurbVarSet<NavierStokes3DVarSet, 0> > >,
		      CellCenterFVMData,
		      FiniteVolumeArcJetModule>
arcJetPhiInsulatedWallTurbFVMCCProvider("ArcJetPhiInsulatedWallTurbFVMCC");

// electrode
MethodCommandProvider<ArcJetPhiElectrode<NoSlipWallIsothermalNSKOmPvt
                      <EulerVarSet, NavierStokesTurbVarSet<NavierStokes3DVarSet, 0> > >,
                      CellCenterFVMData,
                      FiniteVolumeArcJetModule>
arcJetPhiElectrodeTurbFVMCCProvider("ArcJetPhiElectrodeTurbFVMCC");


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
