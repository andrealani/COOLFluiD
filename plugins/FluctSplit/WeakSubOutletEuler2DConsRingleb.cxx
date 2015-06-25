#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubOutletEuler2DConsRingleb.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSubOutletEuler2DConsRingleb, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubOutletEuler2DConsringlebProvider("WeakSubOutletEuler2DConsRingleb");

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsRingleb::defineConfigOptions(Config::OptionList& options)
{
   //options.addConfigOption< CFreal >("alpha","Flow direction");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler2DConsRingleb::WeakSubOutletEuler2DConsRingleb(const std::string& name) :
  WeakBC2D(name),
  _varSet()
{
   addConfigOptionsTo(this);
    //m_alpha = 1.0;
   //setParameter("alpha",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler2DConsRingleb::~WeakSubOutletEuler2DConsRingleb()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsRingleb::setup()
{
  CFAUTOTRACE;

  WeakBC2D::setup();

  _varSet->setup();

 
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsRingleb::setGhostState(const State& state,
					     State& gstate)
{
const CFreal gamma = _varSet->getModel()->getGamma();

  CFreal _angle=90.0;

  
  //First we compute qbar which is the velocity referred to the stagnation speed of sound
 
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal oneovergammaMinus1 = 1./gammaMinus1;
  //First we compute qbar which is the velocity referred to the stagnation speed of sound
 
  CFreal qbar = 0.43;
  CFreal qbar2 = qbar*qbar;
  CFreal cbar = sqrt(1.0 - 0.5*gammaMinus1*qbar2);
  CFreal rhobar =  exp((2.0*oneovergammaMinus1)*log(cbar));
  CFreal rhobar2 = rhobar*rhobar;
  const Node& node = state.getCoordinates();

  //Then we compute the k
  
  CFreal x = node[XX];

  
  CFreal Jbar = (1.0/cbar) + (1.0/(3.0*exp(3.0*log(cbar)))) + (1.0/(5.0*exp(5.0*log(cbar)))) - 0.5*log10((1.0+cbar)/(1.0-cbar));
  CFreal k = sqrt(2.0/((1.0/(qbar2)) - 2.*rhobar*(x - Jbar/2.)));

  CFreal oneoverrho2 = 1.0/rhobar2;
  CFreal oneoverq2 = 1.0/qbar2;
  CFreal oneoverk2 = 1.0/(k*k);
  CFreal oneoverk = 1.0/(k);
  CFreal oneoverrho = 1.0/rhobar;
  CFreal oneoverq = 1.0/qbar;

  CFreal dcbar = -(gammaMinus1*qbar*0.5)/sqrt(1.0-gammaMinus1*0.5*qbar2);

  CFreal drhobar = (2.0*oneovergammaMinus1)*exp((2.0*oneovergammaMinus1-1.0)*log(cbar))*dcbar;

  CFreal dj = -((1.0/(cbar*cbar)) + exp(-4.0*log(cbar)) + exp(-6.0*log(cbar)) +
              (1.0/(log(10.)*(1.0-cbar*cbar))))*dcbar;

  CFreal dx =-( 0.5*(-oneoverrho2*drhobar)*(oneoverq2 - 2.0*oneoverk2) +
              oneoverrho*(-1.0/(qbar2*qbar)) + 0.5*dj);

  CFreal temp = sqrt(1.0-qbar2*oneoverk2);


  CFreal dy = (oneoverk*oneoverq*oneoverrho2*drhobar*temp +
              oneoverk*oneoverrho*oneoverq2*temp + 
              oneoverk*oneoverk2*oneoverrho*(1.0/temp));

  //CF_DEBUG_OBJ(dx);
  ///CF_DEBUG_OBJ(dy);

  CFreal norm = sqrt(dx*dx + dy*dy);

  dx /=norm;
  dy /=norm;

  
  if ((dy >= 0.) && (dx >= 0.))
  _angle = atan(dy/dx);
  
  if ((dy <= 0.) && (dx <= 0.))
  _angle = atan(dy/dx)+MathTools::MathConsts::CFrealPi();

  if ((dy <= 0.) && (dx >= 0.))
  _angle = atan(dy/dx);

  if ((dy >= 0.) && (dx <= 0.))
  _angle = atan(dy/dx)+MathTools::MathConsts::CFrealPi();



  CFreal rho = state[0];
  CFreal u = state[1]/state[0];
  CFreal v = state[2]/state[0];
  CFreal norm_U = sqrt(u*u + v*v);
  
  gstate[0] = rho;
  gstate[1] = rho*norm_U*cos(_angle);
  gstate[2] = rho*norm_U*sin(_angle);
  gstate[3] = state[3];

}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsRingleb::configure ( Config::ConfigArgs& args )
{
  WeakBC2D::configure(args);

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
		 create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
