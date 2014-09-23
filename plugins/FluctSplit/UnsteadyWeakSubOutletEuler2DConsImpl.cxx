#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "UnsteadyWeakSubOutletEuler2DConsImpl.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
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

MethodCommandProvider<UnsteadyWeakSubOutletEuler2DConsImpl, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> unsteadyWeakSubOutletEuler2DConsImplProvider("UnsteadyWeakSubOutletEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DConsImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Function defining the static pressure.");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSubOutletEuler2DConsImpl::UnsteadyWeakSubOutletEuler2DConsImpl(const std::string& name) :
  WeakBC2DImpl(name),
  m_varSet(),
  m_dUoutDu(),
  m_identity()
{
   addConfigOptionsTo(this);
  m_functions = std::vector<std::string>();
   setParameter("Def",&m_functions);

  m_vars = std::vector<std::string>();
   setParameter("Vars",&m_vars);

  m_pressure = 1.0;
  //addConfigOption("P","static pressure",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyWeakSubOutletEuler2DConsImpl::~UnsteadyWeakSubOutletEuler2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DConsImpl::setup()
{
  WeakBC2DImpl::setup();
  m_varSet->setup();

  m_dUoutDu.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_dUoutDu = 0.;
  m_identity.resize(PhysicalModelStack::getActive()->getNbEq());
  m_identity = 1.;

}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DConsImpl::configure ( Config::ConfigArgs& args )
{
  WeakBC2DImpl::configure(args);

  m_vFunction.setFunctions(m_functions);
  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler2DCons";
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(m_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyWeakSubOutletEuler2DConsImpl::computeFluxAndJacob(
                         vector<State*>& states,
             RealVector& flux,
             RealMatrix& fluxJacob)
{
  // the first state is the ghost one
  State *const gstate = states[0];
  State *const state = states[1];

  // evaluate the value of the pressure
  RealVector variables(PhysicalModelStack::getActive()->getDim() + 1);
  const RealVector& temp = state->getCoordinates();
  for (CFuint i = 0; i < temp.size();++i){
    variables[i] = temp[i];
    }
  variables[temp.size()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  RealVector result(1);
  m_vFunction.evaluate(variables,result);
  cf_assert (result.size() == 1);
  m_pressure = result[0]/m_varSet->getModel()->getPressRef();

  // Define the ghost state
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

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
