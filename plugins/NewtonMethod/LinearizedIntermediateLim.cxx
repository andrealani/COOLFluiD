// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "LinearizedIntermediateLim.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/LSSMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LinearizedIntermediateLim, NewtonIteratorData, NewtonMethodModule> LinearizedIntermediateLimProvider("LinearizedIntermediateLim");

//////////////////////////////////////////////////////////////////////////////

void LinearizedIntermediateLim::execute()
{
  // past rhs
  DataHandle<CFreal> pastRhs = socket_pastRhs.getDataHandle();

  // rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  //limiter
  DataHandle<CFreal> limiter = socket_limiter.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = limiter.size()/nbEqs;
  
  // the 0.5 factor CANNOT be put here (it's in the ComputeRHS !!)
  // because there is no way to access the jacobian matrix
  // and multiply for 0.5 also that
  if (  getMethodData().isAchieved() && SubSystemStatusStack::getActive()->isSubIterationLastStep())
  {
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	const CFuint currID = iState*nbEqs + iEq;
        const CFreal lastRhs = rhs[currID];

        rhs[currID] = 2.*rhs[currID] + limiter[currID] * (pastRhs[currID] - rhs[currID]);
	
        // set the past steady RHS (reference RHS for the next
        // physical time step) to the last computed steady RHS
        pastRhs[currID] = lastRhs;
      }
    }
  }
  else {
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	const CFuint currID = iState*nbEqs + iEq;
        rhs[currID] = 2.*rhs[currID] + limiter[currID] * (pastRhs[currID] - rhs[currID]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > LinearizedIntermediateLim::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastRhs);
  result.push_back(&socket_limiter);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
