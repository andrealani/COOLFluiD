#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubInletEuler3DCons.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSubInletEuler3DCons, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubInletEuler3DConsProvider("WeakSubInletEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler3DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("angleZen","angleZen");
   options.addConfigOption< CFreal >("angleAzi","angleAzi");
   options.addConfigOption< CFreal >("Ptot","total pressure");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler3DCons::WeakSubInletEuler3DCons(const std::string& name) :
  WeakBC3D(name),
  m_varSet(),
  m_dUinDu(),
  m_identity()
{
   addConfigOptionsTo(this);
  m_tTotal = 1.0;
   setParameter("Ttot",&m_tTotal);

  m_pTotal = 1.0;
   setParameter("Ptot",&m_pTotal);

  m_angleAzi = 0.0;
   setParameter("angleAzi",&m_angleAzi);
   
  m_angleZen = 0.0;
  setParameter("angleZen", &m_angleZen);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler3DCons::~WeakSubInletEuler3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler3DCons::setup()
{
  WeakBC3D::setup();
  m_varSet->setup();
  m_dUinDu.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_dUinDu = 0.;
  m_identity.resize(PhysicalModelStack::getActive()->getNbEq());
  m_identity = 1.;

  //convert angle from degrees to radiants
  m_angleAzi *= MathTools::MathConsts::CFrealPi()/180.;
  m_angleZen *= MathTools::MathConsts::CFrealPi()/180.;
  m_tTotal /= m_varSet->getModel()->getTempRef();
  m_pTotal /= m_varSet->getModel()->getPressRef();

}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler3DCons::setGhostState(const State& state,
					     State& gstate)
{
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = m_varSet->getModel()->getR();
  const CFreal rho = state[0];
  const CFreal rhoE = state[4];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal w = state[3]/rho;
  const CFreal vel2 = u*u + v*v + w*w;
  const CFreal p = gammaMinus1*(rhoE - 0.5*rho*vel2);
  // mach number is extrapolated from inside the domain
  const CFreal a2 = gamma*p/rho;
  const CFreal mach = sqrt(vel2/a2);
  const CFreal M2 = mach*mach;
  const CFreal coeffM = 1. + 0.5*gammaMinus1*M2;
  const CFreal ghostT = m_tTotal/coeffM;
  const CFreal ghostP = m_pTotal/pow(coeffM, gammaDivGammaMinus1);
  //const CFreal tgAngle = tan(m_angle);
  /*const CFreal tgAngleAzi = tan(m_angleAzi);
  const CFreal tgAngleZen = tan(m_angleZen);
  const CFreal coeffAngleAzi = 1. + tgAngleAzi*tgAngleAzi;
  const CFreal coeffAngleZen = 1. + tgAngleZen*tgAngleZen;*/
  // density
  gstate[0] = ghostP/(R*ghostT);
  // velocity
  gstate[1] = mach*sqrt(gamma*R*ghostT)*sin(m_angleZen)*cos(m_angleAzi)*gstate[0];
  gstate[2] = mach*sqrt(gamma*R*ghostT)*sin(m_angleZen)*sin(m_angleAzi)*gstate[0];
  gstate[3] = mach*sqrt(gamma*R*ghostT)*cos(m_angleZen)*gstate[0];
  gstate[4] = ghostP/gammaMinus1 + 0.5*(gstate[1]*gstate[1] +
                                           gstate[2]*gstate[2] +
                                           gstate[3]*gstate[3])/gstate[0];

  // here you are transforming 2 states but the transformer works
  // on MaxNbStatesInCell states
  
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler3DCons::configure ( Config::ConfigArgs& args )
{
  WeakBC3D::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string varSetName = "Euler3DCons";
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(m_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
