#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCMirrorEuler2D.hh"
#include "FluxReconstructionNavierStokes/NSBoundaryState.hh"

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
    BCMirrorEuler2D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCMirrorEuler2DProvider("MirrorEuler2D");

//////////////////////////////////////////////////////////////////////////////

BCMirrorEuler2D::BCMirrorEuler2D(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorEuler2D::~BCMirrorEuler2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler2D::computeGhostStates(const vector< State* >& intStates,
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
    // normal
    const RealVector& normal = normals[iState];

    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 4);
    cf_assert(ghostState.size() == 4);
    CFLog(VERBOSE, "state: " << intState << "\n");
    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);
    CFLog(VERBOSE, "physData: " << m_intSolPhysData << "\n");

    // compute normal velocity component
    const CFreal uNX2 = 2.0*(m_intSolPhysData[EulerTerm::VX]*normal[XX] +
                             m_intSolPhysData[EulerTerm::VY]*normal[YY]);

    // set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::VX]  = m_intSolPhysData[EulerTerm::VX] - uNX2*normal[XX];
    m_ghostSolPhysData[EulerTerm::VY]  = m_intSolPhysData[EulerTerm::VY] - uNX2*normal[YY];
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

void BCMirrorEuler2D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 4);

    // normal
    const RealVector& normal = normals[iState];

    // tangential unit vector
    RealVector tangent(2);
    tangent[XX] = -normal[YY];
    tangent[YY] =  normal[XX];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][3];
    RealVector& tempGradG = *ghostGrads[iState][3];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];

    // internal normal and tangential component
    const RealVector uNGradI = uGradI*normal [XX] + vGradI*normal [YY];
    const RealVector uTGradI = uGradI*tangent[XX] + vGradI*tangent[YY];

    // ghost normal and tangential component
    const RealVector uNGradG = uNGradI;
    const CFreal nGradUT = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY];
    const RealVector uTGradG = uTGradI - 2.0*nGradUT*normal;

    // project onto x- and y-axis
    uGradG = uNGradG*normal[XX] + uTGradG*tangent[XX];
    vGradG = uNGradG*normal[YY] + uTGradG*tangent[YY];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler2D::computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                         const std::vector< Framework::State* >& intStates,
                                         const std::vector< Framework::State* >& ghostStates,
                                         const std::vector< RealVector >& unitNormals,
                                         const std::vector< RealVector >& flxPntCoords,
                                         std::vector< RealVector* >& bndGradVars)
{
  setSlipWallBndGradVars(gradVarsFlxPnt,unitNormals,m_velocityIDs,bndGradVars);
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler2D::computeBndStates(const std::vector< Framework::State* >& intStates,
                                       const std::vector< Framework::State* >& ghostStates,
                                       const std::vector< RealVector >& unitNormals,
                                       const std::vector< RealVector >& flxPntCoords,
                                       std::vector< RealVector* >& bndStates)
{
  const CFuint nbrStates = intStates.size();

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    m_eulerVarSet->computePhysicalData(*intStates[iState],m_intSolPhysData);

    // interior pressure and temperature
    m_bndPrimState[0] = m_intSolPhysData[EulerTerm::P];
    m_bndPrimState[3] = m_intSolPhysData[EulerTerm::T];

    // tangential velocity u - (u.n) n
    CFreal normalVel = 0.;
    for (CFuint iDim = 0; iDim < 2; ++iDim)
    {
      normalVel += m_intSolPhysData[EulerTerm::VX+iDim]*unitNormals[iState][iDim];
    }

    for (CFuint iDim = 0; iDim < 2; ++iDim)
    {
      m_bndPrimState[1+iDim] = m_intSolPhysData[EulerTerm::VX+iDim] - normalVel*unitNormals[iState][iDim];
    }

    computeNSBoundaryState(*m_eulerVarSet,*intStates[iState],m_bndPrimState,*bndStates[iState]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler2D::computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                                      std::vector< std::vector< RealVector* > >& bndGrads,
                                      const std::vector< RealVector* >& bndStates,
                                      const std::vector< RealVector >& unitNormals,
                                      const std::vector< RealVector >& flxPntCoords)
{
  setSlipWallBndGrads(intGrads,bndGrads,unitNormals,m_velocityIDs);
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in BCMirrorEuler2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // boundary primitive variables
  m_bndPrimState.resize(4);

  // velocity components in the gradient variables
  m_velocityIDs.resize(2);
  m_velocityIDs[0] = 1;
  m_velocityIDs[1] = 2;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

