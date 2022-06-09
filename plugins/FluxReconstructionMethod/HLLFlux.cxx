#include <iterator>

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/BaseTerm.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/HLLFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    HLLFlux,FluxReconstructionSolverData,RiemannFlux,FluxReconstructionModule >
  HLLFluxProvider("HLLFlux");

//////////////////////////////////////////////////////////////////////////////

HLLFlux::HLLFlux(const std::string& name) :
  RiemannFlux(name),
  m_rightFlux(),
  m_leftFlux(),
  m_rightEv(),
  m_leftEv()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

HLLFlux::~HLLFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& HLLFlux::computeFlux(Framework::State& lState,
                                           Framework::State& rState,
                                           const RealVector& normal)
{  
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  //vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();

  // Set members to current left and right update state
  m_updateStates[LEFT]  = &lState;
  m_updateStates[RIGHT] = &rState;

  //CFLog(VERBOSE, "stateLF = "  << rState << "\n");
  // compute physical data for the left and the right internal flux points
  updateVarSet->computePhysicalData(lState, m_pData[LEFT]);
  updateVarSet->computePhysicalData(rState, m_pData[RIGHT]);

  // flux for the right and left state
  //const RealVector& unitNormal = getMethodData().getUnitNormal();
  m_rightFlux = updateVarSet->getFlux()(m_pData[RIGHT], normal);
  m_leftFlux = updateVarSet->getFlux()(m_pData[LEFT], normal);

  // here put update variables
  updateVarSet->computeEigenValues(m_pData[RIGHT], normal, m_rightEv);
  updateVarSet->computeEigenValues(m_pData[LEFT], normal, m_leftEv);

  CFreal aR = 0.0;
  CFreal aL = 0.0;
  for (CFuint i = 0; i < m_nbrEqs; ++i) 
  {
    aR = max(aR, m_rightEv[i]);
    aL = max(aL, m_leftEv[i]);
  }

  const CFreal amax  = max((CFreal)0.0,max(aR,aL));

  aR = 0.0;
  aL = 0.0;
  for (CFuint i = 0; i < m_nbrEqs; ++i) 
  {
    aR = min(aR, m_rightEv[i]);
    aL = min(aL, m_leftEv[i]);
  }

  const CFreal amin  = min((CFreal)0.0,min(aR,aL));
  
  // transform from update states (which are stored) to solution states (in which the equations are written)
  m_solStates[LEFT ] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[LEFT ]); 

  const RealVector lSolState = *(m_solStates)[LEFT ];
  
  m_solStates[RIGHT] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[RIGHT]);

  const RealVector rSolState = *(m_solStates)[RIGHT];

  m_rFlux = 0.5*(m_leftFlux + m_rightFlux -
		(amax + amin)/(amax - amin) * (m_rightFlux - m_leftFlux) +
		amax * amin/(amax - amin) * (rSolState - lSolState));

  return m_rFlux;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& HLLFlux::computeFlux(State& lState, RealVector& lExtraVars,
                                           State& rState, RealVector& rExtraVars,
                                           const RealVector& normal)
{
  // There is no implementation for extravars yet.
  return computeFlux(lState,rState,normal);
}

//////////////////////////////////////////////////////////////////////////////

void HLLFlux::setup()
{
  CFAUTOTRACE;

  RiemannFlux::setup();
  
  m_rightFlux.resize(m_nbrEqs);
  m_leftFlux.resize(m_nbrEqs);
  m_rightEv.resize(m_nbrEqs);
  m_leftEv.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void HLLFlux::unsetup()
{
  CFAUTOTRACE;

  RiemannFlux::unsetup();
  
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
