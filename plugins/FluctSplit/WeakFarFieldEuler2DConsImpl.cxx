#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakFarFieldEuler2DConsImpl.hh"
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

MethodCommandProvider<WeakFarFieldEuler2DConsImpl, FluctuationSplitData, FluctSplitNavierStokesModule> weakFarFieldEuler2DConsImplProvider("WeakFarFieldEuler2DConsImpl");

//////////////////////////////////////////////////////////////////////////////

void WeakFarFieldEuler2DConsImpl::defineConfigOptions(Config::OptionList& options)
{
    options.addConfigOption< CFreal >("Uinf","x velocity");
   options.addConfigOption< CFreal >("Pinf","static pressure");
   options.addConfigOption< CFreal >("Tinf","static temperature");
}

//////////////////////////////////////////////////////////////////////////////

WeakFarFieldEuler2DConsImpl::WeakFarFieldEuler2DConsImpl(const std::string& name) :
  WeakBC2DImpl(name),
  m_varSet(),
  m_dUinDu(),
  m_identity()
{
   addConfigOptionsTo(this);
   
  m_uinf = 0.0;
  setParameter("Uinf",&m_uinf);

  m_pinf = 0.0;
  setParameter("Pinf",&m_pinf);

  m_tinf = 0.0;
  setParameter("Tinf",&m_tinf);




}

//////////////////////////////////////////////////////////////////////////////

WeakFarFieldEuler2DConsImpl::~WeakFarFieldEuler2DConsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakFarFieldEuler2DConsImpl::setup()
{
  WeakBC2DImpl::setup();
  m_varSet->setup();
  m_dUinDu.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_dUinDu = 0.;
  m_identity.resize(PhysicalModelStack::getActive()->getNbEq());
  m_identity = 1.;

  m_pinf /=m_varSet->getModel()->getPressRef();
  m_tinf /=m_varSet->getModel()->getTempRef();
  m_uinf /=m_varSet->getModel()->getVelRef();

}

//////////////////////////////////////////////////////////////////////////////

void WeakFarFieldEuler2DConsImpl::computeFluxAndJacob(
                                     vector<State*>& states,
                                     RealVector& flux,
                                     RealMatrix& fluxJacob)
{  
  State *const gstate = states[0];
  State *const state = states[1];
  // check if we are at inlet or outlet
  CFreal u =  (*state)[1]/(*state)[0];
  CFreal v =  (*state)[2]/(*state)[0];
  const CFreal vel2 = u*u + v*v;

  CFreal ovcoef = 1.0/(sqrt(vel2)*sqrt(m_faceNormal[0]*m_faceNormal[0]+m_faceNormal[1]*m_faceNormal[1]));
  CFreal un =  (u*m_faceNormal[0]+v*m_faceNormal[1])*ovcoef;
  //inlet case un >= 0
    if (un >= 0.0) {

      CFreal cosAngle = u/sqrt(vel2);
      CFreal sinAngle = v/sqrt(vel2);
      const CFreal gamma = m_varSet->getModel()->getGamma();
      const CFreal gammaMinus1 = gamma - 1.;
      const CFreal R = m_varSet->getModel()->getR();
  
      (*gstate)[0] = m_pinf/(R*m_tinf);
      (*gstate)[1] = cosAngle*m_uinf*(*gstate)[0];
      (*gstate)[2] = sinAngle*m_uinf*(*gstate)[0];
      (*gstate)[3] = m_pinf/gammaMinus1 + 0.5*((*gstate)[1]*(*gstate)[1] +
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
  CFreal vel32 = pow((((*state)[1])*((*state)[1]) +
		      ((*state)[2])*((*state)[2])), 3.0/2.0);
  CFreal ovsqrtvel = pow((((*state)[1])*((*state)[1]) +
		      ((*state)[2])*((*state)[2])), -0.5);
  m_dUinDu = 0.;
 
  m_dUinDu(1,1) = m_uinf*(ovsqrtvel - (((*state)[1])*((*state)[1]))/vel32);
  m_dUinDu(1,2) = m_uinf*( -(((*state)[2])*((*state)[1]))/vel32);
  m_dUinDu(1,3) = 0.0;

  m_dUinDu(2,0) = 0.0;
  m_dUinDu(2,1) = m_uinf*( -(((*state)[2])*((*state)[1]))/vel32);
  m_dUinDu(2,2) = m_uinf*(ovsqrtvel - (((*state)[2])*((*state)[2]))/vel32);
  m_dUinDu(2,3) = 0.0;


  m_dUinDu -= m_identity;

  fluxJacob = m_kPlus*m_dUinDu;

   }

  //outlet un <0
        else { 
	  // the first state is the ghost one
	  

	  (*gstate)[0] = (*state)[0];
	  (*gstate)[1] = (*state)[1];
	  (*gstate)[2] = (*state)[2];
	  (*gstate)[3] = m_pinf/(m_varSet->getModel()->getGamma() - 1.) +
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
	  m_dUinDu = 0.;
	  m_dUinDu(3,0) = -0.5*vel2;
	  m_dUinDu(3,1) = u;
	  m_dUinDu(3,2) = v;
	  m_dUinDu(3,3) = -1.;
	  
	  fluxJacob = m_kPlus*m_dUinDu;
	}
}

//////////////////////////////////////////////////////////////////////////////

void WeakFarFieldEuler2DConsImpl::configure ( Config::ConfigArgs& args )
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
