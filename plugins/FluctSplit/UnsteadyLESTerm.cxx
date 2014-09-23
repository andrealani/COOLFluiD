#include "FluctSplit/UnsteadyLESTerm.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"


#include "FluctSplit/FluctSplitLES.hh"
#include "FluctSplit/UnsteadyLESTerm.hh"

#include "LESvki/LESvki.hh"
#include "LESvki/Gradient2DVarSet.hh"
#include "LESvki/Clark2DVarSet.hh"
#include "LESvki/Smagorinsky2DVarSet.hh"
#include "LESvki/WALES2DVarSet.hh"
#include "ComputeDiffusiveTerm.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::LESvki;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluctSplit {

      MethodStrategyProvider<UnsteadyLESTerm<Gradient2DVarSet>,
 			     FluctuationSplitData,
 			     ComputeDiffusiveTerm,
 			     FluctSplitLESModule>
 unsteadygradientDiffusiveTermProvider("UnsteadyGradient");

      MethodStrategyProvider<UnsteadyLESTerm<Smagorinsky2DVarSet>,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitLESModule>
UnsteadysmagDiffusiveTermProvider("UnsteadySmagorinsky");

      MethodStrategyProvider<UnsteadyLESTerm<Clark2DVarSet>,
 		       FluctuationSplitData,
 		       ComputeDiffusiveTerm,
 		       FluctSplitLESModule>
 UnsteadyClarkDiffusiveTermProvider("UnsteadyClark");

  MethodStrategyProvider<UnsteadyLESTerm<WALES2DVarSet>,
	       FluctuationSplitData,
	       ComputeDiffusiveTerm,
	       FluctSplitLESModule>
UnsteadywalesDiffusiveTermProvider("UnsteadyWALES");


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
