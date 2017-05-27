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

LaxFriedrichsFlux::LaxFriedrichsFlux(const std::string& name) :
  RiemannFlux(name),
  m_sumFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedrichsFlux::~LaxFriedrichsFlux()
{
  CFAUTOTRACE;
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

  CFLog(VERBOSE, "stateLF = "  << rState << "\n");
  // compute physical data for the left and the right internal flux points
  updateVarSet->computePhysicalData(lState, m_pData[LEFT]);
  updateVarSet->computePhysicalData(rState, m_pData[RIGHT]);
  CFLog(VERBOSE, "HERE2!!\n");
  
  // flux for right and left state (the physical data must be passed here!)
  m_sumFlux  = updateVarSet->getFlux()(m_pData[LEFT], normal);
  m_sumFlux += updateVarSet->getFlux()(m_pData[RIGHT], normal);
  
  // compute left and right maximum absolute eigenvalues
  const CFreal lMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(m_pData[LEFT], normal);
  const CFreal rMaxAbsEVal = updateVarSet->getMaxAbsEigenValue(m_pData[RIGHT], normal);
  
  // compute absoluteJacobian |A|
  const CFreal absA = 0.5*(lMaxAbsEVal+rMaxAbsEVal);
  
  // transform from update states (which are stored) to solution states (in which the equations are written)
  m_solStates[LEFT ] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[LEFT ]);              
  m_solStates[RIGHT] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[RIGHT]);
  
  State& lSolState = *(m_solStates)[LEFT ];
  State& rSolState = *(m_solStates)[RIGHT];
  
  // compute the Riemann flux
  // Flux = 1/2*(Fmin + Fplus) - 1/2*|A|*(Uplus - Umin)
  m_rFlux = 0.5*(m_sumFlux -  absA*(rSolState - lSolState));
  
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
