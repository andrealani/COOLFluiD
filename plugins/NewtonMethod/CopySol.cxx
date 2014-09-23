// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "CopySol.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
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

MethodCommandProvider<CopySol, NewtonIteratorData, NewtonMethodModule> CopySolProvider("CopySol");

//////////////////////////////////////////////////////////////////////////////

void CopySol::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = states.size();

  for (CFuint iState = 0; iState < nbStates; ++iState) {
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

      CFLogDebugMax( "rhs(iState,iEq) " << rhs(iState, iEq, nbEqs) << "\n");
      (*states[iState])[iEq] = rhs(iState, iEq, nbEqs);

    }
  }

  // reset to 0 the update coefficient
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CopySol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics
} // namespace COOLFluiD
