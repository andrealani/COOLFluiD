#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/BCMirrorEuler3D.hh"

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
    BCMirrorEuler3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionNavierStokesModule >
  BCMirrorEuler3DProvider("MirrorEuler3D");

//////////////////////////////////////////////////////////////////////////////

BCMirrorEuler3D::BCMirrorEuler3D(const std::string& name) :
  BCStateComputer(name),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorEuler3D::~BCMirrorEuler3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler3D::computeGhostStates(const vector< State* >& intStates,
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

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // compute normal velocity component
    const CFreal uNX2 = 2.0*(m_intSolPhysData[EulerTerm::VX]*normal[XX] +
                             m_intSolPhysData[EulerTerm::VY]*normal[YY] +
                             m_intSolPhysData[EulerTerm::VZ]*normal[ZZ]);

    // set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::VX]  = m_intSolPhysData[EulerTerm::VX] - uNX2*normal[XX];
    m_ghostSolPhysData[EulerTerm::VY]  = m_intSolPhysData[EulerTerm::VY] - uNX2*normal[YY];
    m_ghostSolPhysData[EulerTerm::VZ]  = m_intSolPhysData[EulerTerm::VZ] - uNX2*normal[ZZ];
    m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
    m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_intSolPhysData[EulerTerm::V]*
                                                  m_intSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];

    // set the ghost state from its physical data
    m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                            std::vector< std::vector< RealVector* > >& ghostGrads,
                                            const std::vector< RealVector >& normals,
                                            const std::vector< RealVector >& coords)
{
  /// @NOTE this does not work ...

  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 5);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY] + presGradI[ZZ]*normal[ZZ]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][4];
    RealVector& tempGradG = *ghostGrads[iState][4];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY] + tempGradI[ZZ]*normal[ZZ]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];
    RealVector& wGradI = *intGrads  [iState][3];
    RealVector& wGradG = *ghostGrads[iState][3];

    // internal normal and tangential component
    const RealVector velocNGradI = uGradI*normal[XX] + vGradI*normal[YY] + wGradI*normal[ZZ];
    const RealVector uTGradI = uGradI - velocNGradI*normal[XX];
    const RealVector vTGradI = vGradI - velocNGradI*normal[YY];
    const RealVector wTGradI = wGradI - velocNGradI*normal[ZZ];

    // ghost normal and tangential component
    const RealVector velocNGradG = velocNGradI;
    RealVector velocTGradNI(3);
    velocTGradNI[XX] = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY] + uTGradI[ZZ]*normal[ZZ];
    velocTGradNI[YY] = vTGradI[XX]*normal[XX] + vTGradI[YY]*normal[YY] + vTGradI[ZZ]*normal[ZZ];
    velocTGradNI[ZZ] = wTGradI[XX]*normal[XX] + wTGradI[YY]*normal[YY] + wTGradI[ZZ]*normal[ZZ];
    const RealVector uTGradG = uTGradI - 2.0*velocTGradNI[XX]*normal;
    const RealVector vTGradG = vTGradI - 2.0*velocTGradNI[YY]*normal;
    const RealVector wTGradG = wTGradI - 2.0*velocTGradNI[ZZ]*normal;

    // compute ghost velocity gradients
    uGradG = uTGradG + velocNGradG*normal[XX];
    vGradG = vTGradG + velocNGradG*normal[YY];
    wGradG = wTGradG + velocNGradG*normal[ZZ];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorEuler3D::setup()
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
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCMirrorEuler3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

