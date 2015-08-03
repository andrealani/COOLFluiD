#include "FiniteVolumeNavierStokes/FiniteVolumeLTE.hh"
#include "FiniteVolumeNavierStokes/SubBCTurb.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVTLTE.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVT.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "LTE/Euler2DPuvtLTE.hh"
#include "LTE/Euler3DPvtLTE.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LTE;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

// inlet with given mass flow, temperature
MethodCommandProvider<SubBCTurb
		      <SubInletEulerPvtVTLTE, MultiScalarVarSet<Euler2DPuvtLTE>, true>, 
		      CellCenterFVMData, 
		      FiniteVolumeLTEModule> 
subInletTurb2DMassFlowLTEFVMCCProvider("SubInletTurb2DMassFlowLTEFVMCC");

MethodCommandProvider<SubBCTurb
		      <SubInletEulerPvtVTLTE, MultiScalarVarSet<Euler3DPvtLTE>, true>, 
		      CellCenterFVMData, 
		      FiniteVolumeLTEModule> 
subInletTurb3DMassFlowLTEFVMCCProvider("SubInletTurb3DMassFlowLTEFVMCC");

// inlet with given velocity, temperature
MethodCommandProvider<SubBCTurb
		      <SubInletEulerPvtVT<MultiScalarVarSet<Euler2DPuvtLTE> >, 
		       MultiScalarVarSet<Euler2DPuvtLTE>, true>, 
		      CellCenterFVMData, 
		      FiniteVolumeLTEModule> 
subInletTurb2DUVTLTEFVMCCProvider("SubInletTurb2DVTLTEFVMCC");
      
MethodCommandProvider<SubBCTurb
		      <SubInletEulerPvtVT<MultiScalarVarSet<Euler3DPvtLTE> >, 
		       MultiScalarVarSet<Euler3DPvtLTE>, true>, 
		      CellCenterFVMData, 
		      FiniteVolumeLTEModule> 
subInletTurb3DVTLTEFVMCCProvider("SubInletTurb3DVTLTEFVMCC");
      
// outlet with given presssure
MethodCommandProvider<SubBCTurb
		      <SubOutletEulerPvt<MultiScalarVarSet<Euler3DPvtLTE> >, 
		       MultiScalarVarSet<Euler3DPvtLTE>, false>, 
		      CellCenterFVMData, 
		      FiniteVolumeLTEModule> 
subOutletTurb3DVTLTEFVMCCProvider("SubOutletTurb3DVTLTEFVMCC");
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
