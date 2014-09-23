// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "BDF2_CN1stStepPrepare.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BDF2_CN1stStepPrepare, NewtonIteratorData, NewtonMethodModule> BDF2_CN1stStepPrepareProvider("CN1stStepPrepare");

//////////////////////////////////////////////////////////////////////////////

void BDF2_CN1stStepPrepare::execute()
{
  DataHandle<CFreal> rhs         = socket_rhs        .getDataHandle();
  DataHandle<CFreal> pastTimeRhs = socket_pastTimeRhs.getDataHandle();
  cf_assert(rhs.size() == pastTimeRhs.size());
  const CFuint nbrRHS = rhs.size();
  for (CFuint iRHS = 0; iRHS < nbrRHS; ++iRHS)
  {
    pastTimeRhs[iRHS] = rhs[iRHS];
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > BDF2_CN1stStepPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastTimeRhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
