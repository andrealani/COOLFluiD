#include "FiniteVolumeNEQ/ThermNEQST.hh"
#include "FiniteVolumeNEQ/ChemNEQST.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
     
MethodStrategyProvider<ThermNEQST<MultiScalarVarSet<Euler1DVarSet> >,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
euler1DChemThermNEQFVMCCProvider("Euler1DCTNEQST");
 
MethodStrategyProvider<ThermNEQST<MultiScalarVarSet<Euler2DVarSet> >,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
euler2DChemThermNEQFVMCCProvider("Euler2DCTNEQST");

MethodStrategyProvider<ThermNEQST<MultiScalarVarSet<Euler3DVarSet> >,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeNEQModule>
euler3DChemThermNEQFVMCCProvider("Euler3DCTNEQST");

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
