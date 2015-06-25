#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubInletEuler2DConsRingleb.hh"
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

MethodCommandProvider<WeakSubInletEuler2DConsRingleb, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubInletEuler2DConsringlebProvider("WeakSubInletEuler2DConsRingleb");

void WeakSubInletEuler2DConsRingleb::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Mach","Mach number");
   options.addConfigOption< CFreal >("Ptot","total pressure");

}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DConsRingleb::WeakSubInletEuler2DConsRingleb(const std::string& name) :
  WeakBC2D(name),
  _varSet()
{
   addConfigOptionsTo(this);
   m_tTotal = 1.0;
   setParameter("Ttot",&m_tTotal);

  _pTotal = 1.0;
   setParameter("Ptot",&_pTotal);
   m_mach = 1.0;
   setParameter("Mach",&m_mach);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DConsRingleb::~WeakSubInletEuler2DConsRingleb()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsRingleb::setup()
{
  CFAUTOTRACE;

  WeakBC2D::setup();
_varSet->setup();

  _varSet->setup();

  m_tTotal /= _varSet->getModel()->getTempRef();
  _pTotal /= _varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsRingleb::setGhostState(const State& state,
					    State& gstate)
{
 const CFreal gamma = _varSet->getModel()->getGamma();
 const CFreal gammaMinus1 = gamma - 1.;
 const CFreal R = _varSet->getModel()->getR();
  //  const CFreal u = state[1]/state[0];
  //  const CFreal v = state[2]/state[0];
  //  const CFreal vel2 = u*u + v*v;

// const Node& node0 = state.getCoordinates();
//   const CFreal p = gammaMinus1*(state[3] - 0.5*state[0]*vel2);
  // mach number is extrapolated from inside the domain
 // const CFreal mach = sqrt(vel2/(gamma*p/state[0]));
  const CFreal coeffM = 1. + 0.5*gammaMinus1*m_mach*m_mach;
  const CFreal ghostT = m_tTotal/coeffM;
//  const CFreal ghostP = _pTotal/pow(coeffM, gamma/gammaMinus1);
CFreal rho = state[0];
const CFreal ghostP = rho*R*ghostT;

  CFreal angle = getAngle(state);

  const CFreal tgAngle = tan(angle);
// CFuint ID = state.getLocalID();
//if ((ID == 1))
  // {
//CF_DEBUG_OBJ(tgAngle);
//}
  const CFreal coeffAngle = 1. + tgAngle*tgAngle;
//CF_DEBUG_OBJ(1./coeffAngle);
//CF_DEBUG_OBJ(tgAngle);

  

if ((tgAngle <= 0.) && (angle >=0))
{
  // density
  gstate[0] = rho;
  gstate[1] = -m_mach*sqrt(gamma*R*ghostT/coeffAngle)*rho;
  gstate[2] = gstate[1]*tgAngle;
  gstate[3] = ghostP/gammaMinus1 + 0.5*(gstate[1]*gstate[1] +
					gstate[2]*gstate[2])/rho;
}
else
{
  // density
  gstate[0] = rho;
  gstate[1] = m_mach*sqrt(gamma*R*ghostT/coeffAngle)*gstate[0];
  gstate[2] = gstate[1]*tgAngle;
  gstate[3] = ghostP/gammaMinus1 + 0.5*(gstate[1]*gstate[1] +
					gstate[2]*gstate[2])/gstate[0];
}
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsRingleb::configure ( Config::ConfigArgs& args )
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
CFreal WeakSubInletEuler2DConsRingleb::getAngle(const State& state)
{
  CFreal _angle=90.0;
//_angle *= MathTools::MathConsts::CFrealPi()/180.;

// CFuint ID = state.getLocalID();
  //First we compute qbar which is the velocity referred to the stagnation speed of sound
  const CFreal gamma = _varSet->getModel()->getGamma();

  const CFreal gammaMinus1 = gamma - 1.;
 
  const CFreal oneovergammaMinus1 = 1.0/gammaMinus1;
  // mach number is extrapolated from inside the domain
  

//   //First we compute qbar which is the velocity referred to the stagnation speed of sound
//  // CFreal qbar = sqrt(M2/(M2*gammaMinus1*0.5+1.));

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

   CFreal dx = 0.5*(-oneoverrho2*drhobar)*(oneoverq2 - 2.0*oneoverk2) +
               oneoverrho*(-1.0/(qbar2*qbar)) + 0.5*dj;

   CFreal temp = sqrt(1.0-qbar2*oneoverk2);


  CFreal dy = oneoverk*oneoverq*oneoverrho2*drhobar*temp +
               oneoverk*oneoverrho*oneoverq2*temp + 
               oneoverk*oneoverk2*oneoverrho*(1.0/temp);
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
 
  return _angle;
}
//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
