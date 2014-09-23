#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/AddBDF2TimeRHSInGivenCellSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< AddBDF2TimeRHSInGivenCellSpectralFD,
                       SpectralFDMethodData,
                       SpectralFDModule>
AddBDF2TimeRHSInGivenCellSpectralFD("BDF2TimeRHSInGivenCell");

//////////////////////////////////////////////////////////////////////////////

AddBDF2TimeRHSInGivenCellSpectralFD::AddBDF2TimeRHSInGivenCellSpectralFD(const std::string& name) :
  AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD(name),
  socket_pastTimeRhs("pastTimeRhs")
{
}

//////////////////////////////////////////////////////////////////////////////

AddBDF2TimeRHSInGivenCellSpectralFD::~AddBDF2TimeRHSInGivenCellSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void AddBDF2TimeRHSInGivenCellSpectralFD::setup()
{
  AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::setup();
}

//////////////////////////////////////////////////////////////////////////////

void AddBDF2TimeRHSInGivenCellSpectralFD::addTimeResidual()
{
  // get factor for the residual
  const CFreal resFactor = getMethodData().getResFactor();
  const CFreal resFactorP1 = resFactor + 1.;

  // get current cell rhs, past states and past time rhs
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();
  DataHandle< State* > pastStates       = socket_pastStates      .getDataHandle();
  DataHandle< CFreal > pastTimeRhs      = socket_pastTimeRhs     .getDataHandle();

  // add time residual contribution
  const CFuint nbrSolPnts = m_cellStates->size();
  CFuint resIdx = 0;
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get state
    const State& currState = *(*m_cellStates)[iSol];

    // get state ID
    const CFuint stateID = currState.getLocalID();

    // get past state
    const State& pastState = *pastStates[stateID];

    // add contribution to rhs
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resIdx)
    {
      const CFreal dUdt = (currState[iEq] - pastState[iEq])*m_diagValues[iSol];
      rhsCurrStatesSet[resIdx] -= resFactorP1*dUdt - resFactor*pastTimeRhs(stateID,iEq,m_nbrEqs);

      // store dUdt for next time step
      if(SubSystemStatusStack::getActive()->isLastStep())
      {
        pastTimeRhs(stateID,iEq,m_nbrEqs) = dUdt;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
AddBDF2TimeRHSInGivenCellSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::needsSockets();

  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
