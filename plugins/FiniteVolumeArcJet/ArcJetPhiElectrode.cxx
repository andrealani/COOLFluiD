#include "FiniteVolumeArcJet/ArcJetPhiElectrode.hh"
#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVT.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"
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

MethodCommandProvider<ArcJetPhiElectrode<NoSlipWallIsothermalNSvt<MultiScalarTerm<EulerTerm> > >, 
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiElectrodeFVMCCProvider("ArcJetPhiElectrodeFVMCC");

MethodCommandProvider<ArcJetPhiElectrode< SubInletEulerPvtVT<EulerVarSet> >, 
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiElectrodeInletFVMCCProvider("ArcJetPhiElectrodeInletFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
