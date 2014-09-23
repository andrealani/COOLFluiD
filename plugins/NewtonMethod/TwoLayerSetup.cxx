// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"
#include "TwoLayerSetup.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
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

MethodCommandProvider<TwoLayerSetup, NewtonIteratorData, NewtonMethodModule> twoLayerSetupProvider("TwoLayerSetup");

//////////////////////////////////////////////////////////////////////////////

TwoLayerSetup::TwoLayerSetup(std::string name) : 
  StdSetup(name),
  socket_interStates("interStates"),
  socket_interRhs("interRhs"),
  socket_interUpdateCoeff("interUpdateCoeff")  
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerSetup::~TwoLayerSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
TwoLayerSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = 
    StdSetup::providesSockets();
  
  result.push_back(&socket_interStates);
  result.push_back(&socket_interRhs);
  result.push_back(&socket_interUpdateCoeff);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSetup::execute()
{
  CFAUTOTRACE;
  
  // first call to the parent class
  StdSetup::execute();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // States at time step u_n
  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  
  interStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    interStates[i] = new State();
  }
  
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();
  interRhs.resize(nbStates*nbEqs);
  interRhs = 0.0;
  
  const CFuint nbLSS = getMethodData().getLinearSystemSolver().size();
  cf_assert(nbLSS > 0);
  DataHandle<CFreal> interUpdateCoeff = 
    socket_interUpdateCoeff.getDataHandle();
  interUpdateCoeff.resize(nbStates*nbLSS);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
