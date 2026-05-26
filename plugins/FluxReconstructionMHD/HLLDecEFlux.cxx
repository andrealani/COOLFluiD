#include <iterator>

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/BaseTerm.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/HLLDecEFlux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    HLLDecEFlux, FluxReconstructionSolverData, RiemannFlux, FluxReconstructionMHDModule >
  HLLDecEFluxProvider("HLLDecEFlux");

//////////////////////////////////////////////////////////////////////////////

void HLLDecEFlux::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >  ("UseAlpha",        "Use modified HLL diffusion coefficient");
  options.addConfigOption< bool >  ("AddLax",          "Add extra Lax dissipation on the energy slot");
  options.addConfigOption< CFreal >("Extra_kappaLax_E","Magnitude of the extra dissipation on the energy slot");
  options.addConfigOption< bool >  ("2Dornot",         "Use the non-dimensional plasma-beta definition");
}

//////////////////////////////////////////////////////////////////////////////

HLLDecEFlux::HLLDecEFlux(const std::string& name) :
  RiemannFlux(name),
  m_rightFlux(),
  m_leftFlux(),
  m_rightEv(),
  m_leftEv(),
  m_lSolState(),
  m_rSolState()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_useAlpha = false;
  setParameter("UseAlpha", &m_useAlpha);

  m_addLax = false;
  setParameter("AddLax", &m_addLax);

  m_kappaLax_E = 0.125;
  setParameter("Extra_kappaLax_E", &m_kappaLax_E);

  m_2Dornot = false;
  setParameter("2Dornot", &m_2Dornot);
}

//////////////////////////////////////////////////////////////////////////////

HLLDecEFlux::~HLLDecEFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& HLLDecEFlux::computeFlux(Framework::State& lState,
                                     Framework::State& rState,
                                     const RealVector& normal)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();

  m_updateStates[LEFT]  = &lState;
  m_updateStates[RIGHT] = &rState;

  updateVarSet->computePhysicalData(lState, m_pData[LEFT]);
  updateVarSet->computePhysicalData(rState, m_pData[RIGHT]);

  m_rightFlux = updateVarSet->getFlux()(m_pData[RIGHT], normal);
  m_leftFlux  = updateVarSet->getFlux()(m_pData[LEFT],  normal);

  updateVarSet->computeEigenValues(m_pData[RIGHT], normal, m_rightEv);
  updateVarSet->computeEigenValues(m_pData[LEFT],  normal, m_leftEv);

  CFreal aR = 0.0;
  CFreal aL = 0.0;
  for (CFuint i = 0; i < m_nbrEqs; ++i) {
    aR = max(aR, m_rightEv[i]);
    aL = max(aL, m_leftEv[i]);
  }
  const CFreal amax = max((CFreal)0.0, max(aR, aL));

  aR = 0.0;
  aL = 0.0;
  for (CFuint i = 0; i < m_nbrEqs; ++i) {
    aR = min(aR, m_rightEv[i]);
    aL = min(aL, m_leftEv[i]);
  }
  const CFreal amin = min((CFreal)0.0, min(aR, aL));

  // transform update states to solution states (separate buffers required;
  // transform() returns a shared internal buffer overwritten by next call)
  m_solStates[LEFT] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[LEFT]);
  m_lSolState = *(m_solStates)[LEFT];

  m_solStates[RIGHT] = getMethodData().getUpdateToSolutionVecTrans()->transform(m_updateStates[RIGHT]);
  m_rSolState = *(m_solStates)[RIGHT];

  CFreal diffCoef = amax * amin / (amax - amin);
  if (m_useAlpha) {
    const CFreal alpha = m_addLax ? 1.0
                                  : max(std::abs(amin), std::abs(amax)) / (amax - amin);
    diffCoef *= 2.0 * alpha;
  }
  m_rFlux = 0.5*(m_leftFlux + m_rightFlux
                 - (amax + amin)/(amax - amin) * (m_rightFlux - m_leftFlux)
                 + diffCoef * (m_rSolState - m_lSolState));

  if (m_addLax) {
    // assumes MHD3DProjection layout: B at [4,5,6], p at [7]
    const CFreal P0  = 0.03851;    // N/m^2
    const CFreal mu0 = 1.2566e-6;
    const CFreal p_tempL = lState[7];
    const CFreal p_tempR = rState[7];
    const CFreal pthL = p_tempL*P0;
    const CFreal pthR = p_tempR*P0;
    const CFreal BxL = lState[4]; const CFreal ByL = lState[5]; const CFreal BzL = lState[6];
    const CFreal BxR = rState[4]; const CFreal ByR = rState[5]; const CFreal BzR = rState[6];
    const CFreal Bsq2L = BxL*BxL + ByL*ByL + BzL*BzL;
    const CFreal Bsq2R = BxR*BxR + ByR*ByR + BzR*BzR;
    const CFreal pmagL = max(Bsq2L/2.0/mu0, 1.e-16);
    const CFreal pmagR = max(Bsq2R/2.0/mu0, 1.e-16);
    CFreal Q_Betafactor;
    if (m_2Dornot) {
      // 2D path: skip the plasma-beta weighting, apply full extra dissipation.
      Q_Betafactor = 1.0;
    }
    else {
      const CFreal plasmaBetaL = pthL/pmagL;
      const CFreal plasmaBetaR = pthR/pmagR;
      const CFreal fac = 1.0;
      const CFreal Q_BetafactorL = std::tanh(fac/plasmaBetaL/100.0);
      const CFreal Q_BetafactorR = std::tanh(fac/plasmaBetaR/100.0);
      Q_Betafactor = max(Q_BetafactorL, Q_BetafactorR);
    }

    m_rFlux[7] -= m_kappaLax_E*Q_Betafactor*max(amax, std::abs(amin))*(m_rSolState[7] - m_lSolState[7]);
  }

  return m_rFlux;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& HLLDecEFlux::computeFlux(State& lState, RealVector& lExtraVars,
                                     State& rState, RealVector& rExtraVars,
                                     const RealVector& normal)
{
  return computeFlux(lState, rState, normal);
}

//////////////////////////////////////////////////////////////////////////////

void HLLDecEFlux::setup()
{
  CFAUTOTRACE;
  RiemannFlux::setup();

  m_rightFlux.resize(m_nbrEqs);
  m_leftFlux.resize(m_nbrEqs);
  m_rightEv.resize(m_nbrEqs);
  m_leftEv.resize(m_nbrEqs);
  m_lSolState.resize(m_nbrEqs);
  m_rSolState.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
