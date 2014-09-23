// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "CrankNichLimIntermediate.hh"
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

MethodCommandProvider<CrankNichLimIntermediate, NewtonIteratorData, NewtonMethodModule> CrankNichLimIntermediateProvider("CrankNichLimIntermediate");

//////////////////////////////////////////////////////////////////////////////

void CrankNichLimIntermediate::execute()
{
  CFAUTOTRACE;

  // past rhs
  DataHandle<CFreal> pastRhs = socket_pastRhs.getDataHandle();

  // rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // limiter in time
  DataHandle<CFreal> timeLimiter = socket_timeLimiter.getDataHandle();

  const CFuint rhsSize = rhs.size();

  // the 0.5 factor CANNOT be put here (it's in the ComputeRHS !!)
  // because there is no way to access the jacobian matrix
  // and multiply for 0.5 also that
  if (getMethodData().isAchieved() && SubSystemStatusStack::getActive()->isSubIterationLastStep()) {
    for (CFuint i = 0; i < rhsSize; ++i) {
      const CFreal lastRhs = rhs[i];
      rhs[i] += pastRhs[i] * 0.5 * timeLimiter[i];
      // set the past steady RHS (reference RHS for the next
      // physical time step) to the last computed steady RHS
      pastRhs[i] = lastRhs / (1. - 0.5 * timeLimiter[i]);
    }
  }
  else {
    for (CFuint i = 0; i < rhsSize; ++i) {
      rhs[i] += pastRhs[i] * 0.5 * timeLimiter[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > CrankNichLimIntermediate::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = CrankNichIntermediate::needsSockets();

  result.push_back(&socket_timeLimiter);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
