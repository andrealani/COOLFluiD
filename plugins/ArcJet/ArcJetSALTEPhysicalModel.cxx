#include "ArcJet/ArcJetSALTE.hh"
#include "ArcJet/ArcJetLTEPhysicalModel.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

ObjectProvider<ArcJetLTEPhysicalModel<MultiScalarTerm<EulerTerm>, NSTurbTerm, DIM_2D>, 
	       PhysicalModelImpl, ArcJetSALTEModule, 1>
arcJetSALTE2DProvider("ArcJetSALTE2D");

ObjectProvider<ArcJetLTEPhysicalModel<MultiScalarTerm<EulerTerm>, NSTurbTerm, DIM_3D>, 
	       PhysicalModelImpl, ArcJetSALTEModule, 1>
arcJetSALTE3DProvider("ArcJetSALTE3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
