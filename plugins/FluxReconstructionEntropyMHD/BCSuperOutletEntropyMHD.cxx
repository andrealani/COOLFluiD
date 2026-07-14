#include "FluxReconstructionEntropyMHD/FluxReconstructionEntropyMHD.hh"
#include "FluxReconstructionEntropyMHD/BCSuperOutletEntropyMHD.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSuperOutletEntropyMHD, FluxReconstructionSolverData, BCStateComputer, FluxReconstructionEntropyMHDModule>
  BCSuperOutletEntropyMHDProvider("SuperOutletEntropyMHD");

//////////////////////////////////////////////////////////////////////////////

BCSuperOutletEntropyMHD::BCSuperOutletEntropyMHD(const std::string& name)
  : BCStateComputer(name)
{
}

//////////////////////////////////////////////////////////////////////////////

BCSuperOutletEntropyMHD::~BCSuperOutletEntropyMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutletEntropyMHD::setup()
{
  BCStateComputer::setup();

  m_needsSpatCoord = true;
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutletEntropyMHD::computeGhostStates(
    const vector<State*>& intStates,
    vector<State*>& ghostStates,
    const vector<RealVector>& normals,
    const vector<RealVector>& coords)
{
  const CFuint nbrFlxPnts = intStates.size();

  for (CFuint iState = 0; iState < nbrFlxPnts; ++iState)
  {
    State& intState   = *intStates[iState];
    State& ghostState = *ghostStates[iState];

    // Zero-gradient: ghost = interior for all equations except phi
    for (CFuint iEq = 0; iEq < 8; ++iEq)
    {
      ghostState[iEq] = intState[iEq];
    }

    // phi: absorbing BC (sign flip so phi=0 at face)
    ghostState[8] = -intState[8];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutletEntropyMHD::computeGhostGradients(
    const vector<vector<RealVector*> >& intGrads,
    vector<vector<RealVector*> >& ghostGrads,
    const vector<RealVector>& normals,
    const vector<RealVector>& coords)
{
  const CFuint nbrFlxPnts = intGrads.size();

  for (CFuint iState = 0; iState < nbrFlxPnts; ++iState)
  {
    const CFuint nbrVars = intGrads[iState].size();
    for (CFuint iVar = 0; iVar < nbrVars; ++iVar)
    {
      *ghostGrads[iState][iVar] = *intGrads[iState][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
