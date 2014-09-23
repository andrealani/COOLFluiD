#include <numeric>

#include "LinearizedEuler2DSourceTermWall.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplitLinEuler.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler2DSourceTermWall,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler2DSTProviderWall("LinEuler2DSourceWall");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermWall::LinearizedEuler2DSourceTermWall(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  m_alpha(-1.0),
  m_eps(-1.0),
  m_freq(-1.0),
  _sourceloc(0)
{
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermWall::~LinearizedEuler2DSourceTermWall()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermWall::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermWall::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
//IMPORTANT!!!!!
// You should not multiply here by delta t, this is done later for the whole fluctuation!!!!!

  const RealVector centroid = cell->computeCentroid();
  DistributionData& ddata = getMethodData().getDistributionData();
  const CFreal time = ddata.time;
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  vector<State*>& states = *ddata.states;
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal pi=3.14159265;
  const CFreal alpha = 17.32868;
  const CFreal eps = 1.0;
  const CFreal om = 8.0*pi;
  const CFreal xs = 0.0;
  const CFreal ys = 0.0;

  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)  {

    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    const CFreal f=eps*exp(-alpha*((x-xs)*(x-xs)+(y-ys)*(y-ys)));
    
//     cout << "time: " << time;

    source[0] += 0.0;
    source[1] += 0.0;
    source[2] += 0.0;
    source[3] += f*sin(om*time);

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
