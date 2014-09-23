#include "Framework/MethodStrategyProvider.hh"

#include "LinEuler/LinEuler2DVarSet.hh"
#include "LinEuler/LinEulerTerm.hh"

#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"
#include "SpectralFDLinEuler/BCMirrorLinEuler2D.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::LinearizedEuler;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCMirrorLinEuler2D,SpectralFDMethodData,BCStateComputer,SpectralFDLinEulerModule >
  BCMirrorLinEuler2DProvider("MirrorLinEuler2D");

//////////////////////////////////////////////////////////////////////////////

BCMirrorLinEuler2D::BCMirrorLinEuler2D(const std::string& name) :
  BCStateComputer(name),
  m_linEulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorLinEuler2D::~BCMirrorLinEuler2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorLinEuler2D::computeGhostStates(const vector< State* >& intStates,
                                            vector< State* >& ghostStates,
                                            const std::vector< RealVector >& normals,
                                            const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // normal
    const RealVector& normal = normals[iState];

    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size()   == 4);
    cf_assert(ghostState.size() == 4);

    // set the extra variables
    m_linEulerVarSet->setExtraPhysicalVars((*m_extraVars)[iState]);

    // set the physical data starting from the inner state
    m_linEulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // compute normal velocity component
    const CFreal uNX2 = 2.0*(m_intSolPhysData[LinEulerTerm::u]*normal[XX] +
                             m_intSolPhysData[LinEulerTerm::v]*normal[YY]);

    // set the physical data for the ghost state
    m_ghostSolPhysData[LinEulerTerm::rho] = m_intSolPhysData[LinEulerTerm::rho];
    m_ghostSolPhysData[LinEulerTerm::u]   = m_intSolPhysData[LinEulerTerm::u] - uNX2*normal[XX];
    m_ghostSolPhysData[LinEulerTerm::v]   = m_intSolPhysData[LinEulerTerm::v] - uNX2*normal[YY];
    m_ghostSolPhysData[LinEulerTerm::p]   = m_intSolPhysData[LinEulerTerm::p];

    // set the ghost state from its physical data
    m_linEulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorLinEuler2D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                               std::vector< std::vector< RealVector* > >& ghostGrads,
                                               const std::vector< RealVector >& normals,
                                               const std::vector< RealVector >& coords)
{
  throw Common::NotImplementedException(FromHere(),"BCMirrorLinEuler2D::computeGhostGradients(): this should not be needed for linearized Euler!!!");
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorLinEuler2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // extra variables required
  m_needsExtraVars = true;

  // get linearized Euler 2D varset
  m_linEulerVarSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  if (m_linEulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException(FromHere(),"Update variable set is not LinEuler2DVarSet in BCMirrorLinEuler2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_linEulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_linEulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

