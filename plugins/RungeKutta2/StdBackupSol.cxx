// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RungeKutta2/RungeKutta2.hh"
#include "StdBackupSol.hh"
#include "Framework/State.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta2 {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdBackupSol, RK2Data, RungeKutta2Module> stdBackupSolProvider("StdBackupSol");
    
//////////////////////////////////////////////////////////////////////////////
    
void StdBackupSol::execute()
{
  DataHandle< RealVector> u0 = socket_u0.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  for (CFuint i = 0; i < u0.size(); i++) {
    u0[i] = *(states[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdBackupSol::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
 
  result.push_back(&socket_u0);
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta2

  } // namespace Numerics

} // namespace COOLFluiD
