#include <numeric>

#include "LinearizedEuler2DSourceTermDampingMonopole.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplitLinEuler.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/SubSystemStatus.hh"
#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/LinEuler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler2DSourceTermDampingMonopole,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler2DSTProviderDampingMonopole("LinEuler2DSourceDampingMonopole");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermDampingMonopole::LinearizedEuler2DSourceTermDampingMonopole(const std::string& name) :
  ComputeSourceTermFSM(name),
  socket_dampingCoeff("dampingCoeff"),
  socket_volumes("volumes"),
  _varSet(CFNULL),
  _temp(),
  m_alpha(-1.0),
  m_eps(-1.0),
  m_freq(-1.0),
  _sourceloc(0)
{
  addConfigOptionsTo(this);

  m_alpha = log(2.)/5.;
  setParameter("alpha", &m_alpha);

  m_eps = 0.01;
  setParameter("eps", &m_eps);

  m_freq = 30.;
  setParameter("freq", &m_freq);

  setParameter("location", &_sourceloc);
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDampingMonopole::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("alpha","Width of the source.");
  options.addConfigOption< CFreal >("eps","Amplitude of the source.");
  options.addConfigOption< CFreal >("freq","Frequency of the source.");
  options.addConfigOption< std::vector<CFreal> >("location","Location of the source.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedEuler2DSourceTermDampingMonopole::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_dampingCoeff);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermDampingMonopole::~LinearizedEuler2DSourceTermDampingMonopole()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDampingMonopole::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDampingMonopole::computeSourceFSM(Framework::GeometricEntity *const cell,
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
  DataHandle<CFreal> dampingCoeff = socket_dampingCoeff.getDataHandle();

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  const CFreal pi=3.14159265;
  const CFreal alpha = m_alpha;
  const CFreal eps = m_eps;
  const CFreal om = 2.*pi*m_freq;
  const CFreal xs = _sourceloc[0];
  const CFreal ys = _sourceloc[1];

  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)  {
    State& curr_state = *states[iState];
    CFuint IDstate = (states[iState])->getLocalID();

    CFreal& nu = dampingCoeff[IDstate];

    const CFreal rho = curr_state[0];
    const CFreal rho0u = curr_state[1];
    const CFreal rho0v = curr_state[2];
    const CFreal p = curr_state[3];

    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    const CFreal f=eps*exp(-alpha*((x-xs)*(x-xs)+(y-ys)*(y-ys)));

    /// Add contribution of monopole first
    source[0] += f*sin(om*time);
    source[1] += 0.0;
    source[2] += 0.0;
    source[3] += c*c*f*sin(om*time);

    /// Add contribution of damping term afterwards
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
