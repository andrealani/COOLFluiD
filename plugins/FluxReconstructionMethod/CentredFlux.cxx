#include <iterator>

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CFSide.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/CentredFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    CentredFlux,FluxReconstructionSolverData,RiemannFlux,FluxReconstructionModule >
CentredFluxProvider("CentredFlux");

//////////////////////////////////////////////////////////////////////////////

CentredFlux::CentredFlux(const std::string& name) :
  RiemannFlux(name),
  m_sumFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CentredFlux::~CentredFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& CentredFlux::computeFlux(State& lState,
                                     State& rState,
                                     const RealVector& normal)
{
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();
  
  // Set members to current left and right update state
  m_updateStates[LEFT]  = &lState;
  m_updateStates[RIGHT] = &rState;
  CFLog(VERBOSE, "LeftState = " << (*(m_updateStates[RIGHT]->getData()))[0] << "\n");
  
  // compute physical data for the left and the right internal flux points
  updateVarSet->computePhysicalData(lState, m_pData[LEFT]);
  updateVarSet->computePhysicalData(rState, m_pData[RIGHT]);
  
  // flux for right and left state (the physical data must be passed here!)
  m_sumFlux  = updateVarSet->getFlux()(m_pData[LEFT], normal);
  m_sumFlux += updateVarSet->getFlux()(m_pData[RIGHT], normal);
  
  // compute the Centred Riemann flux
  // Flux = 1/2*(Fmin + Fplus)  without diffusive part
  m_rFlux = 0.5*(m_sumFlux);
  
  return m_rFlux;
  
}
    
//////////////////////////////////////////////////////////////////////////////

RealVector& CentredFlux::computeFlux(State& lState, RealVector& lExtraVars,
                                     State& rState, RealVector& rExtraVars,
                                     const RealVector& normal)
{
  // There is no implementation for extravars yet.
  return computeFlux(lState,rState,normal);
}

//////////////////////////////////////////////////////////////////////////////

void CentredFlux::setup()
{
  CFAUTOTRACE;

  RiemannFlux::setup();

  // resize variables
  m_sumFlux.resize(m_nbrEqs);
}

void CentredFlux::unsetup()
{
  CFAUTOTRACE;

  RiemannFlux::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
