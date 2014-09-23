#include "FiniteVolumeArcJet/ArcJetPhiOutlet.hh"
#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallIsothermalNSvt.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
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

MethodCommandProvider<ArcJetPhiOutlet<SubOutletEulerPvt<Euler3DVarSet> >, 
		      CellCenterFVMData, 
		      FiniteVolumeArcJetModule> 
arcJetPhiOutlet3DFVMCCProvider("ArcJetPhiOutlet3DFVMCC");

MethodCommandProvider<ArcJetPhiOutlet<SubOutletEulerPvt<Euler2DVarSet> >,                                                     
                      CellCenterFVMData,
                      FiniteVolumeArcJetModule>
arcJetPhiOutlet2DFVMCCProvider("ArcJetPhiOutlet2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
