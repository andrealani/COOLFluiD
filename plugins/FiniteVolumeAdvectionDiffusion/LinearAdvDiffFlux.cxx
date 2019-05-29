#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolumePoisson/PureDiffFlux.hh"
#include "FiniteVolumeAdvectionDiffusion/FiniteVolumeAdvectionDiffusion.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "FiniteVolume/NSFlux.hh"
#include "LinearAdv/AdvectionDiffusionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::LinearAdv;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSFlux<AdvectionDiffusionVarSet>,
                       CellCenterFVMData,
		       ComputeDiffusiveFlux,
                       FiniteVolumeAdvectionDiffusionModule>
linearAdvDiffusiveFluxProvider("LinearAdvDiff");
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
