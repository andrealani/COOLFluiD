// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RungeKutta/RungeKutta.hh"
#include "StdBackupSol.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdBackupSol, RKData, RungeKuttaModule> stdBackupSolProvider("StdBackupSol");

//////////////////////////////////////////////////////////////////////////////

void StdBackupSol::execute()
{
  DataHandle< RealVector> u0 = socket_u0.getDataHandle();
  DataHandle<RealVector> tempu  = socket_tempu.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const bool isFirstStep = getMethodData().getIsFirstStep();
  const bool isLastStep  = getMethodData().getIsLastStep();

  if (isFirstStep) {
    for (CFuint i = 0; i < u0.size(); i++) {
      u0[i] = *(states[i]);
    }
    for (CFuint i = 0; i < tempu.size(); i++) {
      tempu[i] = *(states[i]);
    }
  }
  
  if (isLastStep)
  {
      for (CFuint i = 0; i < tempu.size(); i++) {
	if (states[i]->isParUpdatable()) {
	  *(states[i]) = tempu[i];
	}
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdBackupSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_u0);
  result.push_back(&socket_tempu);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD
