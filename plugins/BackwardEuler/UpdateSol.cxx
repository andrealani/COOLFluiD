// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/MeshData.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"

#include "BackwardEuler/BackwardEuler.hh"
#include "BackwardEuler/UpdateSol.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace BackwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateSol, BwdEulerData, BackwardEulerModule> updateSolProvider("UpdateSol");

//////////////////////////////////////////////////////////////////////////////

UpdateSol::UpdateSol(const std::string& name) : BwdEulerCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

void UpdateSol::execute()
{
  CFAUTOTRACE;
  
  DataHandle <Framework::State*, Framework::GLOBAL> states = 
    socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  for (CFuint iState = 0; iState < nbStates; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        (*states[iState])[iEq] += dU(iState, iEq, nbEqs);
      }
    }
  }

  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateSol::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace BackwardEuler

  } // namespace Numerics
} // namespace COOLFluiD
