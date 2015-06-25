#include <numeric>
#include <algorithm>

#include "FluctSplitLinEuler.hh"
#include "FluctSplit/LinEuler/WeakSubOutletHedsrtom2DCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "LinEuler/LinearizedEuler.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "LinEuler/LinEuler2DVarSet.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataSocketSink.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider < WeakSubOutletHedstrom2DCons,
                        FluctuationSplitData,
                        FluctSplitLinEulerModule >
theWeakSubOutletHedstrom2DConsProvider("WeakSubOutletHedstrom2DCons");

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletHedstrom2DCons::WeakSubOutletHedstrom2DCons(const std::string& name) :
  WeakBC2D(name),
  _varSet(),
  _bcNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletHedstrom2DCons::~WeakSubOutletHedstrom2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletHedstrom2DCons::setup()
{

  WeakBC2D::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletHedstrom2DCons::setGhostState(const State& state,
					     State& gstate)
{


//   transform the conservative variables to characteristics. The characteristic states are defined through the normal perpendicular to the boundary
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];
  const CFuint ID = state.getLocalID();

  const RealVector* bcNormal = &(*bcNormalsInTrs)[ID];

  const CFreal ncharx = -(*bcNormal)[0];
  const CFreal nchary = -(*bcNormal)[1];

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];
  const CFreal oneoverc = 1./c;

// characteristic variables belonging to the current state
  CFreal entropy=(state[0])-(state[3])*oneoverc*oneoverc;
  CFreal vorticity=(state[1])*nchary-(state[2])*ncharx;
  CFreal aplus=(state[1])*ncharx+(state[2])*nchary + oneoverc*(state[3]);
  CFreal aminus=0.0;

// transforming back
  gstate[0] = -(entropy+0.5*oneoverc*(aplus+aminus));
  gstate[1] = -(nchary*vorticity+0.5*ncharx*(aplus-aminus));
  gstate[2] = -(-ncharx*vorticity+0.5*nchary*(aplus-aminus));
  gstate[3] = -(0.5*c*(aplus+aminus));
/*
  gstate[0] = 0.0;
  gstate[1] = 0.0;
  gstate[2] = 0.0;
  gstate[3] = 0.0;*/

}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletHedstrom2DCons::configure ( Config::ConfigArgs& args )
{
  WeakBC2D::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "LinEuler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler2DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
