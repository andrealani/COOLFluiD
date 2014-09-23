// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "LinearizedIntermediate.hh"
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

MethodCommandProvider<LinearizedIntermediate, NewtonIteratorData, NewtonMethodModule> LinearizedIntermediateProvider("LinearizedIntermediate");

//////////////////////////////////////////////////////////////////////////////

void LinearizedIntermediate::execute()
{
  // past rhs
  DataHandle<CFreal> pastRhs = socket_pastRhs.getDataHandle();

  // rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint rhsSize = rhs.size();

  // the 0.5 factor CANNOT be put here (it's in the ComputeRHS !!)
  // because there is no way to access the jacobian matrix
  // and multiply for 0.5 also that
  if (getMethodData().isAchieved() && SubSystemStatusStack::getActive()->isSubIterationLastStep()) {
    for (CFuint i = 0; i < rhsSize; ++i) {
      const CFreal lastRhs = rhs[i];
      rhs[i] += pastRhs[i];

      // set the past steady RHS (reference RHS for the next
      // physical time step) to the last computed steady RHS
      pastRhs[i] = lastRhs;
    }
  }
  else {
    for (CFuint i = 0; i < rhsSize; ++i) {
      rhs[i] += pastRhs[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > LinearizedIntermediate::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
