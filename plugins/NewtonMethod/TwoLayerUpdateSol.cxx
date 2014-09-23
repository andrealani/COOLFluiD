// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "TwoLayerUpdateSol.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerUpdateSol, NewtonIteratorData, NewtonMethodModule> twoLayerUpdateSolProvider("TwoLayerUpdateSol");

//////////////////////////////////////////////////////////////////////////////

void TwoLayerUpdateSol::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<State*> interStates  = socket_interStates.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();

  // rhs is the temporary placeholder for the dU (intermediate states)
  DataHandle<CFreal>& dU = rhs;

  // interRhs is the temporary placeholder for the dU (intermediate states)
  DataHandle<CFreal>& interdU = interRhs;

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint statesSize = states.size();

  for (CFuint iState = 0; iState < statesSize; ++iState) {
    // do the update only if the state is parallel updatable
    if (states[iState]->isParUpdatable()) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
        (*states[iState])[iEq] += dU(iState, iEq, nbEqs);
        (*interStates[iState])[iEq] += interdU(iState, iEq, nbEqs);
      }
    }
  }

  // reset to 0 the update coefficient
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  DataHandle<CFreal> interUpdateCoeff = socket_interUpdateCoeff.getDataHandle();

  updateCoeff = 0.0;
  interUpdateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > TwoLayerUpdateSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_interStates);
  result.push_back(&socket_rhs);
  result.push_back(&socket_interRhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_interUpdateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
