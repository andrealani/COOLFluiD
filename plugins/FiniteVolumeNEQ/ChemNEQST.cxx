#include "FiniteVolumeNEQ/ChemNEQST.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ChemNEQST<MultiScalarVarSet<Euler1DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
euler1DCNEQSTFVMCCProvider("Euler1DCNEQST");

MethodStrategyProvider<ChemNEQST<MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
euler2DCNEQSTFVMCCProvider("Euler2DCNEQST");

MethodStrategyProvider<ChemNEQST<MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
euler3DCNEQSTFVMCCProvider("Euler3DCNEQST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
