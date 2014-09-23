#include "Framework/MethodCommandProvider.hh"
#include "AeroCoef/AeroCoefFVM.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFluxCC.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFluxCC3D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesSkinFrictionHeatFluxCC3D<NavierStokesSkinFrictionHeatFluxCC>,
		      DataProcessingData,
		      AeroCoefFVMModule>
navierStokesSkinFrictionHeatFluxCC3DProvider("NavierStokesSkinFrictionHeatFluxCC3D");
      
//////////////////////////////////////////////////////////////////////////////

} // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
