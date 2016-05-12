#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolumeNavierStokes/AUSMPlusFlux.hh"
#include "FiniteVolumeNavierStokes/AUSMPlusUpFlux.hh"
#include "FiniteVolume/AUSMFluxALE.hh"

#include "FiniteVolumeCombustion/FiniteVolumeCombustion.hh"

#include "KOmega/EulerKOmegaVarSet.hh"
#include "NEQ/Euler2DNEQRhoivt.hh"
#include "NEQ/Euler3DNEQRhoivt.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/// AUSM+up
MethodStrategyProvider<AUSMPlusUpFlux<EulerKOmegaVarSet<Euler2DNEQRhoivt,2> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusUpCombustion2DProvider("AUSMPlusUpCombustion2D");

MethodStrategyProvider<AUSMPlusUpFlux<EulerKOmegaVarSet<Euler3DNEQRhoivt,2> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusUpCombustion3DProvider("AUSMPlusUpCombustion3D");
  
/// AUSM+
MethodStrategyProvider<AUSMPlusFlux<EulerKOmegaVarSet<Euler2DNEQRhoivt,2> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusCombustion2DProvider("AUSMPlusCombustion2D");
      
MethodStrategyProvider<AUSMPlusFlux<EulerKOmegaVarSet<Euler3DNEQRhoivt,2> >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusCombustion3DProvider("AUSMPlusCombustion3D");

/// AUSM+up ALE 
MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux<EulerKOmegaVarSet<Euler2DNEQRhoivt,2> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusUpCombustion2DALEProvider("AUSMPlusUpCombustion2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusUpFlux<EulerKOmegaVarSet<Euler3DNEQRhoivt,2> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusUpCombustion3DALEProvider("AUSMPlusUpCombustion3DALE");
     
/// AUSM+ ALE 
MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<EulerKOmegaVarSet<Euler2DNEQRhoivt,2> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusCombustion2DALEProvider("AUSMPlusCombustion2DALE");

MethodStrategyProvider<AUSMFluxALE<AUSMPlusFlux<EulerKOmegaVarSet<Euler3DNEQRhoivt,2> > >,
		       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
		       FiniteVolumeCombustionModule>
ausmPlusCombustion3DALEProvider("AUSMPlusCombustion3DALE");
   
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
