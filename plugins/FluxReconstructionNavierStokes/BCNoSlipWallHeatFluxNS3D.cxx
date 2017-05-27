#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCNoSlipWallHeatFluxNS3D.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCNoSlipWallHeatFluxNS3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCNoSlipWallHeatFluxNS3DProvider("NoSlipWallHeatFluxNS3D");

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxNS3D::BCNoSlipWallHeatFluxNS3D(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCNoSlipWallHeatFluxNS3D::~BCNoSlipWallHeatFluxNS3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some physical data from the model
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::VX]  = -m_intSolPhysData[EulerTerm::VX];// negate velocity )-- average = 0
    m_ghostSolPhysData[EulerTerm::VY]  = -m_intSolPhysData[EulerTerm::VY];// negate velocity )-- average = 0
    m_ghostSolPhysData[EulerTerm::VZ]  = -m_intSolPhysData[EulerTerm::VZ];// negate velocity )-- average = 0
    m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
    m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_intSolPhysData[EulerTerm::V]*
                                                  m_intSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::T] = m_intSolPhysData[EulerTerm::T];

    // set the ghost state from its physical data
    m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 5);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
//     RealVector& pGradI = *intGrads  [iState][0];
//     RealVector& pGradG = *ghostGrads[iState][0];
//     const CFreal nPGrad = pGradI[XX]*normal[XX] + pGradI[YY]*normal[YY] + pGradI[ZZ]*normal[ZZ];
//     pGradG = pGradI - 2.0*nPGrad*normal;
    *ghostGrads[iState][0] = *intGrads[iState][0];

    // velocity
    *ghostGrads[iState][1] = *intGrads[iState][1];
    *ghostGrads[iState][2] = *intGrads[iState][2];
    *ghostGrads[iState][3] = *intGrads[iState][3];

    // temperature
    RealVector& tempGradI = *intGrads  [iState][4];
    RealVector& tempGradG = *ghostGrads[iState][4];
    const CFreal nTempGrad = tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY] + tempGradI[ZZ]*normal[ZZ];
    tempGradG = tempGradI - 2.0*nTempGrad*normal; // + heatFlux*normal;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNoSlipWallHeatFluxNS3D::setup()
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
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCNoSlipWallHeatFluxNS3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

