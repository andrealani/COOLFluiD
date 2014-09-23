#include <numeric>

#include "LinearizedEuler2DSourceTermFlap.hh"
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

MethodStrategyProvider<LinearizedEuler2DSourceTermFlap,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler2DSTProviderFlap("LinEuler2DSourceFlap");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermFlap::LinearizedEuler2DSourceTermFlap(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  m_alpha1(-1.0),
  m_eps1(-1.0),
  m_freq1(-1.0),
  m_xs1(0.0),
  m_ys1(0.0),
  m_alpha2(-1.0),
  m_eps2(-1.0),
  m_freq2(-1.0),
  m_xs2(0.0),
  m_ys2(0.0),
  socket_meanflow("meanflow"),
  socket_volumes("volumes")
{
  addConfigOptionsTo(this);

  m_alpha1 = log(2.)/5.;
  setParameter("alpha1", &m_alpha1);

  m_eps1 = 0.01;
  setParameter("eps1", &m_eps1);

  m_freq1 = 60.;
  setParameter("freq1", &m_freq1);

  m_xs1 = 0.;
  setParameter("xs1", &m_xs1);

  m_ys1 = 0.;
  setParameter("ys1", &m_ys1);
  
  m_alpha2 = log(2.)/5.;
  setParameter("alpha2", &m_alpha2);

  m_eps2 = 0.01;
  setParameter("eps2", &m_eps2);

  m_freq2 = 60.;
  setParameter("freq2", &m_freq2);

  m_xs2 = 0.;
  setParameter("xs2", &m_xs2);

  m_ys2 = 0.;
  setParameter("ys2", &m_ys2);
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermFlap::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("alpha1","Width of the source.");
  options.addConfigOption< CFreal >("eps1","Amplitude of the source.");
  options.addConfigOption< CFreal >("freq1","Frequency of the source.");
  options.addConfigOption< CFreal >("xs1","x coordinate of the source.");
  options.addConfigOption< CFreal >("ys1","y coordinate the source.");
  options.addConfigOption< CFreal >("alpha2","Width of the source.");
  options.addConfigOption< CFreal >("eps2","Amplitude of the source.");
  options.addConfigOption< CFreal >("freq2","Frequency of the source.");
  options.addConfigOption< CFreal >("xs2","x coordinate of the source.");
  options.addConfigOption< CFreal >("ys2","y coordinate the source.");  
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermFlap::~LinearizedEuler2DSourceTermFlap()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LinearizedEuler2DSourceTermFlap::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_meanflow);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
void LinearizedEuler2DSourceTermFlap::setup()
{
  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  
}


/////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermFlap::computeSourceFSM(Framework::GeometricEntity *const cell,
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
  const CFreal alpha1 = m_alpha1;
  
  const CFreal eps1 = m_eps1;
  const CFreal om1 = 2.*pi/m_freq1;
  const CFreal xs1 = m_xs1;
  const CFreal ys1 = m_ys1;
  
  const CFreal alpha2 = m_alpha2;
  const CFreal eps2 = m_eps2;
  const CFreal om2 = 2.*pi/m_freq2;
  const CFreal xs2 = m_xs2;
  const CFreal ys2 = m_ys2;  

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


    source[0] += 0.0;
    if(r>0.004)
	{
	source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+0.0;
	}
    else
	{
	const CFreal f1=eps1*cos(0.1*pi*(x-xs1))*exp(-alpha1*((y-ys1)*(y-ys1))*1000000.);
 	const CFreal f2=eps2*cos(0.1*pi*(x-xs2))*exp(-alpha2*((y-ys2)*(y-ys2))*1000000.);   
	source[1] += (rho*u0+rho0u)*du0dx + (rho*v0+rho0v)*du0dy+f1*sin(om1*time)+f2*sin(om2*time);
	}
    source[2] += (rho*u0+rho0u)*dv0dx + (rho*v0+rho0v)*dv0dy+0.0;
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
