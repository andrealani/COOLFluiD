#include <iterator>

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/BaseTerm.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/LaxFriedrichsFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    LaxFriedrichsFlux,FluxReconstructionSolverData,RiemannFlux,FluxReconstructionModule >
  LaxFriedrichsFluxProvider("LaxFriedrichsFlux");

//////////////////////////////////////////////////////////////////////////////

void LaxFriedrichsFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Epsilon","Diffusion reduction coefficient for tuning the dissipative part of the Rusanov flux (default = 0.5)");
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedrichsFlux::LaxFriedrichsFlux(const std::string& name) :
  RiemannFlux(name),
  m_sumFlux(),
  m_epsilon(0.5)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  setParameter("Epsilon",&m_epsilon);
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedrichsFlux::~LaxFriedrichsFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedrichsFlux::configure ( Config::ConfigArgs& args )
{
  RiemannFlux::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

RealVector& LaxFriedrichsFlux::computeFlux(Framework::State& lState,
                                           Framework::State& rState,
                                           const RealVector& normal)
{
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();

  // Set members to current left and right update state
  m_updateStates[LEFT]  = &lState;
  m_updateStates[RIGHT] = &rState;

  //CFLog(VERBOSE, "stateLF = "  << rState << "\n");
  // compute physical data for the left and the right internal flux points
  updateVarSet->computePhysicalData(lState, m_pData[LEFT]);
  updateVarSet->computePhysicalData(rState, m_pData[RIGHT]);
  
  // flux for right and left state (the physical data must be passed here!)
  m_sumFlux  = updateVarSet->getFlux()(m_pData[LEFT], normal);
  m_sumFlux += updateVarSet->getFlux()(m_pData[RIGHT], normal);
  
  // compute left and right maximum absolute eigenvalues
  const CFreal lMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(m_pData[LEFT], normal);
  const CFreal rMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(m_pData[RIGHT], normal);
  
  // compute absoluteJacobian |A|
  const CFreal absA = 0.5*(lMaxAbsEVal+rMaxAbsEVal);
  
  // transform from update states (which are stored) to solution states (in which the equations are written)
  // Note: transform() returns a pointer to a shared internal buffer that is overwritten on the next call,
  // so we must copy the left result into pre-allocated storage before computing the right transform.
  m_solStates[LEFT ] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[LEFT ]);
  m_lSolState = *(m_solStates)[LEFT ];

  m_solStates[RIGHT] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[RIGHT]);
  m_rSolState = *(m_solStates)[RIGHT];

  // compute the Riemann flux (Rusanov flux)
  // Flux = 1/2*(Fmin + Fplus) - epsilon*|A|*(Uplus - Umin)
  // where epsilon is the diffusion reduction coefficient
  m_rFlux = 0.5*m_sumFlux - m_epsilon*absA*(m_rSolState - m_lSolState);

  return m_rFlux;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& LaxFriedrichsFlux::computeFlux(State& lState, RealVector& lExtraVars,
                                           State& rState, RealVector& rExtraVars,
                                           const RealVector& normal)
{
  // There is no implementation for extravars yet.
  return computeFlux(lState,rState,normal);
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedrichsFlux::setup()
{
  CFAUTOTRACE;

  RiemannFlux::setup();

  // resize variables
  m_sumFlux.resize(m_nbrEqs);
  m_lSolState.resize(m_nbrEqs);
  m_rSolState.resize(m_nbrEqs);

}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedrichsFlux::unsetup()
{
  CFAUTOTRACE;

  RiemannFlux::unsetup();
  
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
