// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "UFEMPrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
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

MethodCommandProvider<UFEMPrepare, NewtonIteratorData, NewtonMethodModule> UFEMPrepareProvider ( "UFEMPrepare" );

//////////////////////////////////////////////////////////////////////////////

UFEMPrepare::UFEMPrepare ( std::string name ) : NewtonIteratorCom ( name ),
  socket_interStates ( "interStates" ),
  socket_states ( "states" ),
  socket_pastStates ( "pastStates" ),
  socket_pastpastStates ( "pastpastStates" )
{
}

//////////////////////////////////////////////////////////////////////////////

UFEMPrepare::~UFEMPrepare()
{
}

//////////////////////////////////////////////////////////////////////////////

void UFEMPrepare::execute()
{
  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<State*> pastpastStates = socket_pastpastStates.getDataHandle();

  const CFuint nbStates = states.size();
  // const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  CFint OExtrap = 2;
  CFreal c=0.5; //this is ExtrapForward
  //CFreal Dt = 0.1;
  CFreal c2 = 0.;
  CFreal c1 = 0.;
  CFreal c0 = 0.;

  if ( OExtrap == 0 )
  {
    c0 = 1.;
  }
  if ( OExtrap == 1 )
  {
    c0 = c+1.;
    c1 = -c;
  }
  if ( OExtrap == 2 )
  {
    c0 = 1.5*c+.5*c*c+1.;
    c1 = -2.*c-c*c;
    c2 = .5*c+.5*c*c;
  }

  // Set Initial States to current states
  for ( CFuint i = 0; i < nbStates; ++i )
  {
//    *(interStates[i]) = c0*(*(states[i])) +c1*(*(pastStates[i])) +c2*(*(pastpastStates[i]));
    *(interStates[i]) = 1.5*(*(states[i])) -0.5*(*(pastStates[i]));
    *(pastpastStates[i]) = *(pastStates[i]);
    *(pastStates[i]) = *(states[i]);
  }

}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > UFEMPrepare::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back ( &socket_interStates );
  result.push_back ( &socket_states );
  result.push_back ( &socket_pastStates );
  result.push_back ( &socket_pastpastStates );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
