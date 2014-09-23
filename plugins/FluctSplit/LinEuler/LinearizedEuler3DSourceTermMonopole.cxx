#include <numeric>

#include "LinearizedEuler3DSourceTermMonopole.hh"
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

MethodStrategyProvider<LinearizedEuler3DSourceTermMonopole,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler3DSTProviderMonopole("LinEuler3DSourceMonopole");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler3DSourceTermMonopole::LinearizedEuler3DSourceTermMonopole(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  m_alpha(-1.0),
  m_eps(-1.0),
  m_freq(-1.0)
{
  addConfigOptionsTo(this);

  m_alpha = log(2.)/5.;
  setParameter("alpha", &m_alpha);

  m_eps = 0.01;
  setParameter("eps", &m_eps);

  m_freq = 30.;
  setParameter("freq", &m_freq);

  _sourceloc = vector<CFreal>();
   setParameter("location", &_sourceloc);
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceTermMonopole::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("alpha","Width of the source.");
  options.addConfigOption< CFreal >("eps","Amplitude of the source.");
  options.addConfigOption< CFreal >("freq","Frequency of the source.");
  options.addConfigOption< std::vector<CFreal> >("location","Location of the source.");
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler3DSourceTermMonopole::~LinearizedEuler3DSourceTermMonopole()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceTermMonopole::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler3DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler3DSourceTermMonopole::computeSourceFSM(Framework::GeometricEntity *const cell,
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
  const CFreal alpha = m_alpha;
  const CFreal eps = m_eps;
  const CFreal om = 2.*pi/m_freq;
  const CFreal xs = _sourceloc[0];
  const CFreal ys = _sourceloc[1];
  const CFreal zs = _sourceloc[2];

// CF_DEBUG_OBJ(_sourceloc.size());

  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;
  source[4] = 0.0;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)  {

    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    const CFreal z = (states[iState]->getCoordinates())[ZZ];
    const CFreal f=eps*exp(-alpha*((x-xs)*(x-xs)+(y-ys)*(y-ys)+(z-zs)*(z-zs)));

    source[0] += f*sin(om*time);
    source[1] += 0.0;
    source[2] += 0.0;
    source[3] += 0.0;
    source[4] += c*c*f*sin(om*time);

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
