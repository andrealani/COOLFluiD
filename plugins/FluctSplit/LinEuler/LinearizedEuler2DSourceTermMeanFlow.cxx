#include <numeric>

#include "LinearizedEuler2DSourceTermMeanFlow.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplitLinEuler.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/LinEuler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler2DSourceTermMeanFlow,
           FluctuationSplitData,
           ComputeSourceTermFSM,
           FluctSplitLinEulerModule>
linEuler2DSTProviderMeanFlow("LinEuler2DSourceMeanFlow");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermMeanFlow::LinearizedEuler2DSourceTermMeanFlow(const std::string& name) :
  ComputeSourceTermFSM(name),
  socket_meanflow("meanflow"),
  socket_volumes("volumes"),
  _varSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedEuler2DSourceTermMeanFlow::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_meanflow);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermMeanFlow::~LinearizedEuler2DSourceTermMeanFlow()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMeanFlow::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMeanFlow::computeSourceFSM(Framework::GeometricEntity *const cell,
RealVector& source,
const FluctSplit::InwardNormalsData& normalsData)

{
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  vector<State*>& states = *ddata.states;

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal gamma = linearData[LinEulerTerm::GAMMA];


  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;

  CFreal du0dx = 0.0;
  CFreal du0dy = 0.0;
  CFreal dv0dx = 0.0;
  CFreal dv0dy = 0.0;
  CFreal dp0dx = 0.0;
  CFreal dp0dy = 0.0;


  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    CFuint IDstate = (states[iState])->getLocalID();
    RealVector& meanflow_state = meanflow[IDstate];

    const CFreal u0 = meanflow_state[1];
    const CFreal v0 = meanflow_state[2];
    const CFreal p0 = meanflow_state[3];

    du0dx += normalsData.getNodalNormComp(iState,0)*u0;
    du0dy += normalsData.getNodalNormComp(iState,1)*u0;
    dv0dx += normalsData.getNodalNormComp(iState,0)*v0;
    dv0dy += normalsData.getNodalNormComp(iState,1)*v0;
    dp0dx += normalsData.getNodalNormComp(iState,0)*p0;
    dp0dy += normalsData.getNodalNormComp(iState,1)*p0;
  }

   CFreal oneover2vol = 1./(2.0*volumes[cellID]);
   du0dx *= oneover2vol;
   du0dy *= oneover2vol;
   dv0dx *= oneover2vol;
   dv0dy *= oneover2vol;
   dp0dx *= oneover2vol;
   dp0dy *= oneover2vol;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)  {
    State& curr_state = *states[iState];
    CFuint IDstate = (states[iState])->getLocalID();
    RealVector& meanflow_state = meanflow[IDstate];
 
    const CFreal rho0 = meanflow_state[0];
    const CFreal u0 = meanflow_state[1];
    const CFreal v0 = meanflow_state[2];
    
    const CFreal rho = curr_state[0];
    const CFreal rho0u = curr_state[1];
    const CFreal rho0v = curr_state[2];
    const CFreal p = curr_state[3];

    source[0] += 0.0;
    source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy;
    source[2] += (rho*u0+rho0u)*dv0dx + (rho*v0+rho0v)*dv0dy;
    source[3] += (gamma-1.0)*p*(du0dx+dv0dy)-(gamma-1.0)*(rho0u/rho0*dp0dx+rho0v/rho0*dp0dy);
}

  const CFuint dimSource = source.size();

  for (CFuint iSource = 0; iSource < dimSource; ++iSource) {
    source[iSource] /= nbStatesInCell;
    source[iSource] *= volumes[cellID]*dt;
  }


}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
