#include <numeric>

#include "LinearizedEuler2DSourceTermDampingZone.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplitLinEuler.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/LinEuler2DVarSet.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
//using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler2DSourceTermDampingZone,
           FluctuationSplitData,
           ComputeSourceTermFSM,
           FluctSplitLinEulerModule>
linEuler2DSTProviderDampingZone("LinEuler2DSourceDampingZone");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermDampingZone::LinearizedEuler2DSourceTermDampingZone(const std::string& name) :
  ComputeSourceTermFSM(name),
  socket_dampingCoeff("dampingCoeff"),
  socket_volumes("volumes"),
  _varSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedEuler2DSourceTermDampingZone::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_dampingCoeff);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermDampingZone::~LinearizedEuler2DSourceTermDampingZone()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDampingZone::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDampingZone::computeSourceFSM(Framework::GeometricEntity *const cell,
RealVector& source,
const FluctSplit::InwardNormalsData& normalsData)

{
  const RealVector centroid = cell->computeCentroid();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  vector<State*>& states = *ddata.states;

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> dampingCoeff = socket_dampingCoeff.getDataHandle();

  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    State& curr_state = *states[iState];
    CFuint IDstate = (states[iState])->getLocalID();

    CFreal& nu = dampingCoeff[IDstate];

    const CFreal rho = curr_state[0];
    const CFreal rho0u = curr_state[1];
    const CFreal rho0v = curr_state[2];
    const CFreal p = curr_state[3];

    source[0] += -nu * (rho - 0.0) ;
    source[1] += -nu * (rho0u - 0.0);
    source[2] += -nu * (rho0v -0.0);
    source[3] += -nu * (p -0.0);
}

  const CFuint dimSource = source.size();

  for (CFuint iSource = 0; iSource < dimSource; ++iSource) {
    source[iSource] /= nbStatesInCell;
    source[iSource] *= volumes[cellID]*dt;
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace FluctSplit

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
