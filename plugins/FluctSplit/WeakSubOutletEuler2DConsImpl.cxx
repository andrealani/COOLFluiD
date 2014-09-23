#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubOutletEuler2DConsImpl.hh"
#include "InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
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

MethodCommandProvider<WeakSubOutletEuler2DConsImpl, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubOutletEuler2DConsImplProvider("WeakSubOutletEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler2DConsImpl::WeakSubOutletEuler2DConsImpl(const std::string& name) :
  WeakBC2DImpl(name),
  m_varSet(),
  m_dUoutDu(),
  m_identity()
{
   addConfigOptionsTo(this);
  m_pressure = 1.0;
   setParameter("P",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler2DConsImpl::~WeakSubOutletEuler2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsImpl::setup()
{
  CFAUTOTRACE;

  WeakBC2DImpl::setup();

  m_varSet->setup();
  m_dUoutDu.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_dUoutDu = 0.;
  m_identity.resize(PhysicalModelStack::getActive()->getNbEq());
  m_identity = 1.;

  m_pressure /= m_varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsImpl::computeFluxAndJacob(
                         vector<State*>& states,
             RealVector& flux,
             RealMatrix& fluxJacob)
{
  // the first state is the ghost one
  State *const gstate = states[0];
  State *const state = states[1];

  (*gstate)[0] = (*state)[0];
  (*gstate)[1] = (*state)[1];
  (*gstate)[2] = (*state)[2];
  (*gstate)[3] = m_pressure/(m_varSet->getModel()->getGamma() - 1.) +
    0.5*((*gstate)[1]*(*gstate)[1] + (*gstate)[2]*(*gstate)[2])/(*gstate)[0];

   // here you are transforming 2 states but the transformer works
  // on MaxNbStatesInCell states
  vector<State*> *const linearStates =
    m_updateToLinearVecTrans->transform(&states);

  m_twoStates[0] = (*linearStates)[0];
  m_twoStates[1] = (*linearStates)[1];

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
  flux = m_kPlus*(*gstate - *state);

  const CFreal rho = (*state)[0];
  const CFreal u = (*state)[1]/rho;
  const CFreal v = (*state)[2]/rho;
  const CFreal vel2 = u*u + v*v;

  m_dUoutDu = 0.;

  m_dUoutDu(3,0) = -0.5*vel2;
  m_dUoutDu(3,1) = u;
  m_dUoutDu(3,2) = v;
  m_dUoutDu(3,3) = -1.;

  fluxJacob = m_kPlus*m_dUoutDu;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler2DConsImpl::configure ( Config::ConfigArgs& args )
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

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
