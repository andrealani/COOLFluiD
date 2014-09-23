#include "ArcJet/ArcJet.hh"
#include "ArcJetPhysicalModel.hh"
#include "EulerArcJetPhysicalModel.hh"
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

// full model 2D
Environment::ObjectProvider<ArcJetPhysicalModel<DIM_2D>, PhysicalModelImpl, ArcJetModule, 1>
arcJet2DProvider("ArcJet2D");

// full model 3D
Environment::ObjectProvider<ArcJetPhysicalModel<DIM_3D>, PhysicalModelImpl, ArcJetModule, 1>
arcJet3DProvider("ArcJet3D");

// inviscid model 2D
Environment::ObjectProvider<EulerArcJetPhysicalModel<DIM_2D>, PhysicalModelImpl, ArcJetModule, 1>
eulerArcJet2DProvider("EulerArcJet2D");

// inviscid model 3D
Environment::ObjectProvider<EulerArcJetPhysicalModel<DIM_3D>, PhysicalModelImpl, ArcJetModule, 1>
eulerArcJet3DProvider("EulerArcJet3D");

//////////////////////////////////////////////////////////////////////////////

    } // namespace ArcJet

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
