#include "FiniteVolume/AUSMFluxALE.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusFlux.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux.hh"
#include "FiniteVolumeNavierStokes/AUSMLiouSteffenFlux.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
/// AUSM+up

MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux<Euler1DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUp1dALEProvider("AUSMPlusUp1DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux<Euler2DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUp2dALEProvider("AUSMPlusUp2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux<Euler3DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUp3dALEProvider("AUSMPlusUp3DALE");
      
MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux <MultiScalarVarSet<Euler1DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMS1dALEProvider("AUSMPlusUpMS1DALE");


MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux <MultiScalarVarSet<Euler2DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMS2dALEProvider("AUSMPlusUpMS2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux <MultiScalarVarSet<Euler3DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMS3dALEProvider("AUSMPlusUpMS3DALE");

/// AUSM+
MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<Euler1DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlus1dALEProvider("AUSMPlus1DALE");
      
MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<Euler2DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlus2dALEProvider("AUSMPlus2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<Euler3DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlus3dALEProvider("AUSMPlus3DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<MultiScalarVarSet<Euler1DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusMS1dALEProvider("AUSMPlusMS1DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<MultiScalarVarSet<Euler2DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusMS2dALEProvider("AUSMPlusMS2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<MultiScalarVarSet<Euler3DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusMS3dALEProvider("AUSMPlusMS3DALE");

/// AUSM Liou-Steffen

MethodStrategyProvider<AUSMFluxALE<AUSMLiouSteffenFlux<Euler1DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffen1dALEProvider("AUSMLiouSteffen1DALE");
      
MethodStrategyProvider<AUSMFluxALE<AUSMLiouSteffenFlux<Euler2DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffen2dALEProvider("AUSMLiouSteffen2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMLiouSteffenFlux<Euler3DVarSet> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffen3dALEProvider("AUSMLiouSteffen3DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMLiouSteffenFlux
		       <MultiScalarVarSet<Euler1DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffenMS1dALEProvider("AUSMLiouSteffenMS1DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMLiouSteffenFlux<MultiScalarVarSet<Euler2DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffenMS2dALEProvider("AUSMLiouSteffenMS2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMLiouSteffenFlux<MultiScalarVarSet<Euler3DVarSet> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmLiouSteffenMS3dALEProvider("AUSMLiouSteffenMS3DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux<MultiScalarVarSet<Euler2DPuvt> > >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpTurb2DALEProvider("AUSMPlusUpTurb2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux<MultiScalarVarSet<Euler3DPvt<Euler3DVarSet> > > >,
		       CellCenterFVMData,
		       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeNavierStokesModule>
ausmPlusUpMUS3DALEProvider("AUSMPlusUpTurb3DALE");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
