#include "Framework/MethodCommandProvider.hh"
#include "FluctSplit/FluctuationSplit.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "AeroCoef/AeroCoefFS.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFluxFS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::FluctSplit;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesSkinFrictionHeatFluxFS<Euler2DVarSet>,
		      DataProcessingData,
		      AeroCoefFSModule>
navierStokesSkinFrictionHeatFluxFSProvider
("NavierStokesSkinFrictionHeatFluxFS");

MethodCommandProvider<NavierStokesSkinFrictionHeatFluxFS<Euler3DVarSet>,
		      DataProcessingData,
		      AeroCoefFSModule>
navierStokes3DSkinFrictionHeatFluxFSProvider
("NavierStokes3DSkinFrictionHeatFluxFS");



MethodCommandProvider<NavierStokesSkinFrictionHeatFluxFS
		      <MultiScalarVarSet<Euler2DVarSet> >,
		      DataProcessingData,
                      AeroCoefFSModule>
navierStokesMSSkinFrictionHeatFluxFSProvider
("NavierStokesSkinFrictionHeatFluxFSMS");
  
//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////



