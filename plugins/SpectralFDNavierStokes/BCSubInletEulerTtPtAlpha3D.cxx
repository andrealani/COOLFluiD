#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/BCSubInletEulerTtPtAlpha3D.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSubInletEulerTtPtAlpha3D,SpectralFDMethodData,BCStateComputer,SpectralFDNavierStokesModule >
  BCSubInletEulerTtPtAlpha3DProvider("SubInletEulerTtPtAlpha3D");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerTtPtAlpha3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Ttot","total temperature");
  options.addConfigOption< CFreal >("Ptot","total pressure");
  options.addConfigOption< CFreal >("alphaXY","alphaXY");
  options.addConfigOption< CFreal >("alphaXZ","alphaXZ");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletEulerTtPtAlpha3D::BCSubInletEulerTtPtAlpha3D(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_tTotal(),
  m_pTotal(),
  m_alphaXY(),
  m_alphaXZ()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_tTotal = 0.0;
   setParameter("Ttot",&m_tTotal);

  m_pTotal = 0.0;
   setParameter("Ptot",&m_pTotal);

  m_alphaXY = 0.0;
   setParameter("alphaXY",&m_alphaXY);

  m_alphaXZ = 0.0;
   setParameter("alphaXZ",&m_alphaXZ);
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletEulerTtPtAlpha3D::~BCSubInletEulerTtPtAlpha3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerTtPtAlpha3D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some data from the physical model
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma -1.0;
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFreal idGassConst = m_eulerVarSet->getModel()->getR();

  const CFreal tgAlphaXY = tan(m_alphaXY);
  const CFreal tgAlphaXZ = tan(m_alphaXZ);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // compute some helper variables
    const CFreal machInt     = m_intSolPhysData[EulerTerm::V]/m_intSolPhysData[EulerTerm::A];
    const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInt*machInt;
    const CFreal coeffPow    = pow(coefficient, gammaDivGammaMinus1);
    const CFreal tTotalInt   = m_intSolPhysData[EulerTerm::T]*coefficient;
    const CFreal pTotalInt   = m_intSolPhysData[EulerTerm::P]*coeffPow;
    const CFreal tgAlphaXYInt  = m_intSolPhysData[EulerTerm::VY]/m_intSolPhysData[EulerTerm::VX];
    const CFreal tgAlphaXZInt  = m_intSolPhysData[EulerTerm::VZ]/m_intSolPhysData[EulerTerm::VX];

    // set ghost state quantities
    const CFreal tTotalGhost  = 2.0*m_tTotal   - tTotalInt;
    const CFreal pTotalGhost  = 2.0*m_pTotal   - pTotalInt;
    const CFreal tgAlphaXYGhost = 2.0*tgAlphaXY - tgAlphaXYInt;
    const CFreal tgAlphaXZGhost = 2.0*tgAlphaXZ - tgAlphaXZInt;
    const CFreal machGhost    = machInt;
    const CFreal tGhost       = tTotalGhost/coefficient;
    const CFreal pGhost       = pTotalGhost/coeffPow;

    //set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::P]   = pGhost;
    m_ghostSolPhysData[EulerTerm::RHO] = pGhost/(idGassConst*tGhost);
    m_ghostSolPhysData[EulerTerm::VX]  =
        machGhost*sqrt(gamma*idGassConst*tGhost/(1.0 + tgAlphaXYGhost*tgAlphaXYGhost + tgAlphaXZGhost*tgAlphaXZGhost));
    m_ghostSolPhysData[EulerTerm::VY]  = tgAlphaXYGhost*m_ghostSolPhysData[EulerTerm::VX];
    m_ghostSolPhysData[EulerTerm::VZ]  = tgAlphaXZGhost*m_ghostSolPhysData[EulerTerm::VX];
    m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*pGhost
                                       + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                         (m_ghostSolPhysData[EulerTerm::VX]*m_ghostSolPhysData[EulerTerm::VX] +
                                          m_ghostSolPhysData[EulerTerm::VY]*m_ghostSolPhysData[EulerTerm::VY] +
                                          m_ghostSolPhysData[EulerTerm::VZ]*m_ghostSolPhysData[EulerTerm::VZ])
                                         )/m_ghostSolPhysData[EulerTerm::RHO];

    // set the ghost state from its physical data
    m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerTtPtAlpha3D::computeGhostGradients
                                                    (const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  const CFuint nbrGradVars = intGrads[0].size();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletEulerTtPtAlpha3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 3D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCSubInletEulerTtPtAlpha3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // non-dimensionalize pressure and temperature
  m_tTotal /= m_eulerVarSet->getModel()->getTempRef ();
  m_pTotal /= m_eulerVarSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
