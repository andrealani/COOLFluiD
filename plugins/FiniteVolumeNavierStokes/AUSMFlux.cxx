#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolumeNavierStokes/AUSMPlusFlux.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux.hh"
#include "FiniteVolumeNavierStokes/AUSMLowMlimit.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusFlux_Mp.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusFlux_MpW.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux_Mp.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux_MpW.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusUpIcpFlux_cp.hh"
#include "FiniteVolumeNavierStokes/AUSMLiouSteffenFlux.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/Euler3DPvt.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

// AUSM+up
MethodStrategyProvider<AUSMPlusUpFlux<Euler1DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUp1dProvider("AUSMPlusUp1D");

MethodStrategyProvider<AUSMPlusUpFlux<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUp2dProvider("AUSMPlusUp2D");

MethodStrategyProvider<AUSMPlusUpFlux<Euler3DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUp3dProvider("AUSMPlusUp3D");

// AUSM+up multi-scalar
MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler1DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMS1dProvider("AUSMPlusUpMS1D");

MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMS2dProvider("AUSMPlusUpMS2D");
    
MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMS3dProvider("AUSMPlusUpMS3D");

MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler2DPuvt> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMSPvt2dProvider("AUSMPlusUpMSPvt2D");

MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler3DPvt<Euler3DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMSPvt3dProvider("AUSMPlusUpMSPvt3D");

// kept for backward compatibility  
MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler2DPuvt> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpTurb2dProvider("AUSMPlusUpTurb2D");

// kept for backward compatibility  
MethodStrategyProvider<AUSMPlusUpFlux <MultiScalarVarSet<Euler3DPvt<Euler3DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpTurb3dProvider("AUSMPlusUpTurb3D");

// VDH's variants 
MethodStrategyProvider<AUSMLowMlimit<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLowMlimitProvider("AUSMLowMlimit");

MethodStrategyProvider<AUSMPlusFlux_Mp<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusFlux_MpProvider("AUSMPlus_Mp");

MethodStrategyProvider<AUSMPlusFlux_MpW<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusFlux_MpWProvider("AUSMPlus_MpW");

MethodStrategyProvider<AUSMPlusUpFlux_Mp<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpFlux_MpProvider("AUSMPlusUp_Mp");

MethodStrategyProvider<AUSMPlusUpFlux_MpW<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpFlux_MpWProvider("AUSMPlusUp_MpW");

// AUSMPlusUpIcp_cp (only 2D for now)
MethodStrategyProvider<AUSMPlusUpIcpFlux_cp<Euler2DVarSet>,	
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpIcpFlux_cpProvider("AUSMPlusUpIcp_cp");

MethodStrategyProvider<AUSMPlusUpIcpFlux_cp<MultiScalarVarSet<Euler2DVarSet> >,	
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpIcpFlux_cpMS2DProvider("AUSMPlusUpIcp_cpMS2D");

/// AUSM+
MethodStrategyProvider<AUSMPlusFlux<Euler1DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlus1dProvider("AUSMPlus1D");
      
MethodStrategyProvider<AUSMPlusFlux<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlus2dProvider("AUSMPlus2D");

MethodStrategyProvider<AUSMPlusFlux<Euler3DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlus3dProvider("AUSMPlus3D");

MethodStrategyProvider<AUSMPlusFlux<MultiScalarVarSet<Euler1DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusMS1dProvider("AUSMPlusMS1D");

MethodStrategyProvider<AUSMPlusFlux<MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusMS2dProvider("AUSMPlusMS2D");

MethodStrategyProvider<AUSMPlusFlux<MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusMS3dProvider("AUSMPlusMS3D");

/// AUSM Liou-Steffen

MethodStrategyProvider<AUSMLiouSteffenFlux<Euler1DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffen1dProvider("AUSMLiouSteffen1D");
      
MethodStrategyProvider<AUSMLiouSteffenFlux<Euler2DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffen2dProvider("AUSMLiouSteffen2D");

MethodStrategyProvider<AUSMLiouSteffenFlux<Euler3DVarSet>,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffen3dProvider("AUSMLiouSteffen3D");

MethodStrategyProvider<AUSMLiouSteffenFlux
		       <MultiScalarVarSet<Euler1DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffenMS1dProvider("AUSMLiouSteffenMS1D");

MethodStrategyProvider<AUSMLiouSteffenFlux<MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffenMS2dProvider("AUSMLiouSteffenMS2D");

MethodStrategyProvider<AUSMLiouSteffenFlux<MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffenMS3dProvider("AUSMLiouSteffenMS3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
