#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSubOutletEuler3DConsImpl.hh"
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

MethodCommandProvider<WeakSubOutletEuler3DConsImpl, FluctuationSplitData, FluctSplitNavierStokesModule> weakSubOutletEuler3DConsImplProvider("WeakSubOutletEuler3DConsImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler3DConsImpl::WeakSubOutletEuler3DConsImpl(const std::string& name) :
  WeakBC3DImpl(name),
  m_varSet(),
  m_dUoutDu(),
  m_identity()
{
   addConfigOptionsTo(this);
  m_pressure = 1.0;
   setParameter("P",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

WeakSubOutletEuler3DConsImpl::~WeakSubOutletEuler3DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsImpl::setup()
{
  CFAUTOTRACE;

  WeakBC3DImpl::setup();

  m_varSet->setup();
  m_dUoutDu.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_dUoutDu = 0.;
  m_identity.resize(PhysicalModelStack::getActive()->getNbEq());
  m_identity = 1.;

  m_pressure /= m_varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsImpl::computeFluxAndJacob
(vector<State*>& states,
 RealVector& flux,
 RealMatrix& fluxJacob)
{
  // the first state is the ghost one
  State *const gstate = states[0];
  State *const state = states[1];

  (*gstate)[0] = (*state)[0];
  (*gstate)[1] = (*state)[1];
  (*gstate)[2] = (*state)[2];
  (*gstate)[3] = (*state)[3];
  (*gstate)[4] = m_pressure/(m_varSet->getModel()->getGamma() - 1.) +
    0.5*((*gstate)[1]*(*gstate)[1] + (*gstate)[2]*(*gstate)[2] +
	 (*gstate)[3]*(*gstate)[3])/(*gstate)[0];

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
  const CFreal w = (*state)[3]/rho;
  const CFreal vel2 = u*u + v*v + w*w;

  m_dUoutDu = 0.;

  m_dUoutDu(4,0) = -0.5*vel2;
  m_dUoutDu(4,1) = u;
  m_dUoutDu(4,2) = v;
  m_dUoutDu(4,3) = w;
  m_dUoutDu(4,4) = -1.;

  fluxJacob = m_kPlus*m_dUoutDu;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSubOutletEuler3DConsImpl::configure ( Config::ConfigArgs& args )
{
  WeakBC3DImpl::configure(args);

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const std::string varSetName = "Euler3DCons";
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(m_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
