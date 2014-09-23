#include <numeric>

#include "LinearizedEuler2DSourceTermMeanMonopole.hh"
#include "FluctSplit/InwardNormalsData.hh"
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
using namespace COOLFluiD::Physics::LinearizedEuler;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler2DSourceTermMeanMonopole,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler2DSTProviderMeanMonopole("LinEuler2DSourceMeanMonopole");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermMeanMonopole::LinearizedEuler2DSourceTermMeanMonopole(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  socket_meanflow("meanflow"),
  socket_volumes("volumes")
{
  addConfigOptionsTo(this);

  m_alpha = log(2.)/5.;
  setParameter("alpha", &m_alpha);

  m_eps = 0.01;
  setParameter("eps", &m_eps);

  m_freq = 30.;
  setParameter("freq", &m_freq);

  m_xs = 0.;
  setParameter("xs", &m_xs);

  m_ys = 0.;
  setParameter("ys", &m_ys);

}
//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMeanMonopole::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("alpha","Width of the source.");
  options.addConfigOption< CFreal >("eps","Amplitude of the source.");
  options.addConfigOption< CFreal >("freq","Frequency of the source.");
  options.addConfigOption< CFreal >("xs","x coordinate of the source.");
  options.addConfigOption< CFreal >("ys","y coordinate the source.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedEuler2DSourceTermMeanMonopole::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_meanflow);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermMeanMonopole::~LinearizedEuler2DSourceTermMeanMonopole()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMeanMonopole::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMeanMonopole::computeSourceFSM(Framework::GeometricEntity *const cell,
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
  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal gamma = linearData[LinEulerTerm::GAMMA];
  const CFreal c     = linearData[LinEulerTerm::c];


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

  CFreal du0dx = 0.0;
  CFreal du0dy = 0.0;
  CFreal dv0dx = 0.0;
  CFreal dv0dy = 0.0;
  CFreal dp0dx = 0.0;
  CFreal dp0dy = 0.0;

/* compute derivatives of the mean flow */

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

    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];

    const CFreal r = x*x+y*y;

    const CFreal rho0 = meanflow_state[0];
    const CFreal u0 = meanflow_state[1];
    const CFreal v0 = meanflow_state[2];

    const CFreal rho = curr_state[0];
    const CFreal rho0u = curr_state[1];
    const CFreal rho0v = curr_state[2];
    const CFreal p = curr_state[3];

/* monopole */
    const CFreal f=eps*exp(-alpha*((x-xs)*(x-xs)+(y-ys)*(y-ys)));
    source[0] += 0.0+f*sin(om*time);
    source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+0.0;
    source[2] += (rho*u0+rho0u)*dv0dx + (rho*v0+rho0v)*dv0dy+0.0;
    source[3] += (gamma-1.0)*p*(du0dx+dv0dy)-(gamma-1.0)*(rho0u/rho0*dp0dx+rho0v/rho0*dp0dy)+c*c*f*sin(om*time);

/* dipole */
/*    source[0] += 0.0;
    if(r>25.)
	{
	source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+0.0;
	}
    else
	{
	const CFreal f=eps*cos(0.1*pi*x)*exp(-alpha*(y*y));
	source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+f*sin(om*time);
	}
    source[2] += (rho*u0+rho0u)*dv0dx + (rho*v0+rho0v)*dv0dy+0.0;
    source[3] += (gamma-1.0)*p*(du0dx+dv0dy)-(gamma-1.0)*(rho0u/rho0*dp0dx+rho0v/rho0*dp0dy);
*/
// /* quadrupole */
//     source[0] += 0.0;
//     if(r>100.)
// 	{
// 	source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+0.0;
// 	source[2] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+0.0;
// 	}
//     else
// 	{
// 	const CFreal f1=eps*sin(0.05*pi*x)*exp(-alpha*(y*y));
// 	const CFreal f2=-eps*sin(0.05*pi*y)*exp(-alpha*(x*x));
// 	source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+f1*sin(om*time);
// 	source[2] += (rho*u0+rho0u)*dv0dx + (rho*v0+rho0v)*dv0dy+f2*sin(om*time);
// 	}
//     source[3] += (gamma-1.0)*p*(du0dx+dv0dy)-(gamma-1.0)*(rho0u/rho0*dp0dx+rho0v/rho0*dp0dy)+0.0;

 }

  const CFuint dimSource = source.size();

  for (CFuint iSource = 0; iSource < dimSource; ++iSource) {
    source[iSource] /= nbStatesInCell;
    source[iSource] *= volumes[cellID];
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
