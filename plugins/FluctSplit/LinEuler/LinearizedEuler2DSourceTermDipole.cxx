#include <numeric>

#include "LinearizedEuler2DSourceTermDipole.hh"
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

MethodStrategyProvider<LinearizedEuler2DSourceTermDipole,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler2DSTProviderDipole("LinEuler2DSourceDipole");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermDipole::LinearizedEuler2DSourceTermDipole(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  m_alpha(-1.0),
  m_eps(-1.0),
  m_freq(-1.0),
  m_xs(0.0),
  m_ys(0.0)
{
  addConfigOptionsTo(this);

  m_alpha = log(2.)/5.;
  setParameter("alpha", &m_alpha);

  m_eps = 0.01;
  setParameter("eps", &m_eps);

  m_freq = 60.;
  setParameter("freq", &m_freq);

  m_xs = 0.;
  setParameter("xs", &m_xs);

  m_ys = 0.;
  setParameter("ys", &m_ys);
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDipole::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("alpha","Width of the source.");
  options.addConfigOption< CFreal >("eps","Amplitude of the source.");
  options.addConfigOption< CFreal >("freq","Frequency of the source.");
  options.addConfigOption< CFreal >("xs","x coordinate of the source.");
  options.addConfigOption< CFreal >("ys","y coordinate the source.");
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermDipole::~LinearizedEuler2DSourceTermDipole()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDipole::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  
}


/////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermDipole::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
//IMPORTANT!!!!!
// You should not multiply here by delta t, this is done later for the whole fluctuation!!!!!

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFreal time = ddata.time;
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  vector<State*>& states = *ddata.states;
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();


  const CFreal pi=3.14159265;
  const CFreal alpha = m_alpha;
  const CFreal eps = m_eps;
  const CFreal om = 2.*pi/m_freq;
  const CFreal xs = m_xs;
  const CFreal ys = m_ys;

  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;



  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)  {

    const CFreal xreal = (states[iState]->getCoordinates())[XX];
    const CFreal yreal = (states[iState]->getCoordinates())[YY];
    const CFreal r = sqrt((xreal-xs)*(xreal-xs)+(yreal-ys)*(yreal-ys));
    const CFreal f=eps*cos(0.1*pi*(xreal-xs))*exp(-alpha*((yreal-ys)*(yreal-ys)));

    source[0] += 0.0;
    if(r>5.)
	{
	source[1] += 0.0;
	}
    else
	{
	source[1] +=f*sin(om*time);
	}
    source[2] += 0.0;
    source[3] += 0.0;
  }


  const CFuint dimSource = source.size();

  for (CFuint iSource = 0; iSource < dimSource; ++iSource) {
    source[iSource] /= nbStatesInCell;
    source[iSource] *= volumes[cellID]*dt;
    //*dt;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
