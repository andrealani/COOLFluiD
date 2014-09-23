#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubInletEuler2DConsImpl.hh"
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

MethodCommandProvider<WeakSubInletEuler2DConsImpl, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubInletEuler2DConsImplProvider("WeakSubInletEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("angle","angle");
   options.addConfigOption< CFreal >("Ptot","total pressure");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DConsImpl::WeakSubInletEuler2DConsImpl(const std::string& name) :
  WeakBC2DImpl(name),
  m_varSet(),
  m_dUinDu(),
  m_identity()
{
   addConfigOptionsTo(this);
  m_tTotal = 1.0;
   setParameter("Ttot",&m_tTotal);

  m_pTotal = 1.0;
   setParameter("Ptot",&m_pTotal);

  m_angle = 0.0;
   setParameter("angle",&m_angle);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubInletEuler2DConsImpl::~WeakSubInletEuler2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsImpl::setup()
{
  WeakBC2DImpl::setup();
  m_varSet->setup();
  m_dUinDu.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_dUinDu = 0.;
  m_identity.resize(PhysicalModelStack::getActive()->getNbEq());
  m_identity = 1.;

  //convert angle from degrees to radiants
  m_angle *= MathTools::MathConsts::CFrealPi()/180.;
  m_tTotal /= m_varSet->getModel()->getTempRef();
  m_pTotal /= m_varSet->getModel()->getPressRef();

}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsImpl::computeFluxAndJacob(
                                     vector<State*>& states,
                                     RealVector& flux,
                                     RealMatrix& fluxJacob)
{
  // the first state is the ghost one
  State *const gstate = states[0];
  State *const state = states[1];
  
  const CFreal angle = getAngle(state[1]);
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = m_varSet->getModel()->getR();
  const CFreal rho = (*state)[0];
  const CFreal rhoE = (*state)[3];
  const CFreal u = (*state)[1]/rho;
  const CFreal v = (*state)[2]/rho;
  const CFreal vel2 = u*u + v*v;
  const CFreal p = gammaMinus1*(rhoE - 0.5*rho*vel2);
  // mach number is extrapolated from inside the domain
  const CFreal a2 = gamma*p/rho;
  const CFreal mach = sqrt(vel2/a2);
  const CFreal M2 = mach*mach;
  const CFreal coeffM = 1. + 0.5*gammaMinus1*M2;
  const CFreal ghostT = m_tTotal/coeffM;
  const CFreal ghostP = m_pTotal/pow(coeffM, gammaDivGammaMinus1);
  const CFreal tgAngle = tan(angle);
  const CFreal coeffAngle = 1. + tgAngle*tgAngle;
  // density
  (*gstate)[0] = ghostP/(R*ghostT);
  (*gstate)[1] = mach*sqrt(gamma*R*ghostT/coeffAngle)*(*gstate)[0];
  (*gstate)[2] = (*gstate)[1]*tgAngle;
  (*gstate)[3] = ghostP/gammaMinus1 + 0.5*((*gstate)[1]*(*gstate)[1] +
                                           (*gstate)[2]*(*gstate)[2])/(*gstate)[0];

  // here you are transforming 2 states but the transformer works
  // on MaxNbStatesInCell states
  vector<State*> *const linearStates =
    m_updateToLinearVecTrans->transform(&states);

  m_twoStates[0] = (*linearStates)[0];
  m_twoStates[1] = (*linearStates)[1];

  // linearize the states in the cell
  m_linearizer->linearize(m_twoStates);

  const CFreal kCoeff = 1./ PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // adimensionalize the normal
  const CFreal faceArea = m_faceNormal.norm2();
  m_adimNormal = m_faceNormal/faceArea;

  // compute the Kplus in conservative variables
  m_varSet->computeEigenValuesVectors(m_rightEv,
                                m_leftEv,
                                m_eValues,
                                m_adimNormal);

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    m_eValuesP[iEq] = max(0.,m_eValues[iEq]);
  }
  m_kPlus = m_rightEv*(m_eValuesP*m_leftEv);
  m_kPlus *= kCoeff*faceArea;

  // flux in conservative variables
  flux = m_kPlus*((*gstate) - (*state));

  const CFreal ghostV = sqrt((*gstate)[1]*(*gstate)[1] +
                             (*gstate)[2]*(*gstate)[2])/(*gstate)[0];

  const CFreal abv1 = rho*a2*coeffM;
  const CFreal abv2 = rho*(2.*a2 + gammaMinus1*vel2);
  const CFreal abv3 = ghostV*abv2;
  const CFreal abv4 = gamma*gamma*gammaMinus1*ghostP;
  const CFreal abv5 = gamma*gammaMinus1*(*gstate)[0];
  const CFreal abv6 = gamma*gamma*ghostP;
  const CFreal cosAlpha = cos(angle);
  const CFreal sinAlpha = sin(angle);

  const CFreal drhovDrho  = -abv4*rhoE*M2*(1.-M2)/(rho*abv3);
  const CFreal drhovDrhou = 2.*abv4*u*rhoE*(1.-M2)/(rho*a2*abv3);
  const CFreal drhovDrhov = 2.*abv4*v*rhoE*(1.-M2)/(rho*a2*abv3);
  const CFreal drhovDrhoE = -abv4*M2*(1.-M2)/abv3;

  m_dUinDu(0,0) = 0.5*abv5*rhoE*M2/(abv1*rho);
  m_dUinDu(0,1) = -abv5*rhoE*u/(abv1*rho*a2);
  m_dUinDu(0,2) = -abv5*rhoE*v/(abv1*rho*a2);
  m_dUinDu(0,3) = 0.5*abv5*M2/abv1;

  m_dUinDu(1,0) = drhovDrho*cosAlpha;
  m_dUinDu(1,1) = drhovDrhou*cosAlpha;
  m_dUinDu(1,2) = drhovDrhov*cosAlpha;
  m_dUinDu(1,3) = drhovDrhoE*cosAlpha;

  m_dUinDu(2,0) = drhovDrho*sinAlpha;
  m_dUinDu(2,1) = drhovDrhou*sinAlpha;
  m_dUinDu(2,2) = drhovDrhov*sinAlpha;
  m_dUinDu(2,3) = drhovDrhoE*sinAlpha;

  m_dUinDu(3,0) = abv6*rhoE*M2*(1.-gammaMinus1*(1.-0.5*M2))/(rho*abv2);
  m_dUinDu(3,1) = abv6*rhoE*u*(gammaMinus1*(2.-M2)-2.)/(rho*a2*abv2);
  m_dUinDu(3,2) = abv6*rhoE*v*(gammaMinus1*(2.-M2)-2.)/(rho*a2*abv2);
  m_dUinDu(3,3) = abv6*M2*(1.-gammaMinus1*(1.-0.5*M2))/abv2;

  m_dUinDu -= m_identity;

  fluxJacob = m_kPlus*m_dUinDu;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubInletEuler2DConsImpl::configure ( Config::ConfigArgs& args )
{
  WeakBC2DImpl::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler2DCons";
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(m_varSet.isNotNull());

}
//////////////////////////////////////////////////////////////////////////////
CFreal WeakSubInletEuler2DConsImpl::getAngle(const State& state)
{
  return m_angle;
}
//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
