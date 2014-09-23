#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/AddBDF3TimeRHSInGivenCellSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< AddBDF3TimeRHSInGivenCellSpectralFD,
                       SpectralFDMethodData,
                       SpectralFDModule>
AddBDF3TimeRHSInGivenCellSpectralFD("BDF3TimeRHSInGivenCell");

//////////////////////////////////////////////////////////////////////////////

AddBDF3TimeRHSInGivenCellSpectralFD::AddBDF3TimeRHSInGivenCellSpectralFD(const std::string& name) :
  AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD(name),
  socket_pastTimeRhs("pastTimeRhs"),
  socket_pastPastTimeRhs("pastPastTimeRhs"),
  m_3StepsTMSparams()
{
}

//////////////////////////////////////////////////////////////////////////////

AddBDF3TimeRHSInGivenCellSpectralFD::~AddBDF3TimeRHSInGivenCellSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void AddBDF3TimeRHSInGivenCellSpectralFD::setup()
{
  AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::setup();

  // get the 3 steps time marching scheme parameters (BDF3 with variable time step)
  m_3StepsTMSparams = getMethodData().get3StepsTMSParams();
}

//////////////////////////////////////////////////////////////////////////////

void AddBDF3TimeRHSInGivenCellSpectralFD::addTimeResidual()
{
  // get current cell rhs, past states, past past states and current residual in cell
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();
  DataHandle< State* > pastStates       = socket_pastStates      .getDataHandle();
  DataHandle< CFreal > pastTimeRhs      = socket_pastTimeRhs     .getDataHandle();
  DataHandle< CFreal > pastPastTimeRhs  = socket_pastPastTimeRhs .getDataHandle();

  // Derefence and resize the 3 steps time marching scheme parameters (BDF3 with variable time step)
  RealVector& params3StepsTMS = *m_3StepsTMSparams;

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
      rhsCurrStatesSet[resIdx] -= params3StepsTMS[0]*dUdt
                                - params3StepsTMS[1]*m_diagValues[iSol]*pastState[iEq]
                                - params3StepsTMS[2]*(m_diagValues[iSol]*pastState[iEq] - pastTimeRhs(stateID,iEq,m_nbrEqs))
                                - params3StepsTMS[3]*(m_diagValues[iSol]*pastState[iEq] - pastTimeRhs(stateID,iEq,m_nbrEqs) - pastPastTimeRhs(stateID,iEq,m_nbrEqs));

      // store dUdt at the two previous time steps, for next time step
      if(SubSystemStatusStack::getActive()->isLastStep())
      {
        pastTimeRhs(stateID,iEq,m_nbrEqs) = dUdt;
        pastPastTimeRhs(stateID,iEq,m_nbrEqs) = pastTimeRhs(stateID,iEq,m_nbrEqs);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
AddBDF3TimeRHSInGivenCellSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::needsSockets();

  result.push_back(&socket_pastTimeRhs);
  result.push_back(&socket_pastPastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
