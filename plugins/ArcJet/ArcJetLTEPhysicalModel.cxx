#include "ArcJet/ArcJet.hh"
#include "ArcJetLTEPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ArcJet {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ArcJetLTEPhysicalModel<DIM_2D>, PhysicalModelImpl, ArcJetModule, 1>
arcJetLTE2DProvider("ArcJetLTE2D");

Environment::ObjectProvider<ArcJetLTEPhysicalModel<DIM_3D>, PhysicalModelImpl, ArcJetModule, 1>
arcJetLTE3DProvider("ArcJetLTE3D");

Environment::ObjectProvider<ArcJetLTEPhysicalModel<DIM_2D>, PhysicalModelImpl, ArcJetModule, 1>
arcJetPG2DProvider("ArcJetPG2D");

Environment::ObjectProvider<ArcJetLTEPhysicalModel<DIM_3D>, PhysicalModelImpl, ArcJetModule, 1>
arcJetPG3DProvider("ArcJetPG3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
