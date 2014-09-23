// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ForwardEuler/ForwardEuler.hh"
#include "TwoLayerSetup.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerSetup, FwdEulerData, ForwardEulerLib> twoLayerSetupProvider("TwoLayerSetup");

TwoLayerSetup::TwoLayerSetup(std::string name) :
  FwdEulerCom(name),
  socket_pastStates("pastStates"),
  socket_interStates("interStates"),
  socket_rhs("rhs"),
  socket_interRhs("interRhs"),
  socket_updateCoeff("updateCoeff"),
  socket_interUpdateCoeff("interUpdateCoeff"),
  socket_states("states")
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
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);
  result.push_back(&socket_rhs);
  result.push_back(&socket_interRhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_interUpdateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  getMethodData().setNormSize(2);

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  if (SubSystemStatusStack::getActive()->getDT() > 0.){
    // States at time step u_n
    DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

    if(pastStates.size() != nbStates)
    {
      pastStates.resize(nbStates);
      for (CFuint i = 0; i < nbStates; ++i) {
        pastStates[i] = new State();
      }
    }
  }

  // States at time step u_n
  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  interStates.resize(nbStates);
  for (CFuint i = 0; i < nbStates; ++i) {
    interStates[i] = new State();
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(nbStates*nbEqs);
  rhs = 0.0;

  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();
  interRhs.resize(nbStates*nbEqs);
  interRhs = 0.0;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(nbStates);

  DataHandle<CFreal> interUpdateCoeff = socket_interUpdateCoeff.getDataHandle();
  interUpdateCoeff.resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > TwoLayerSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD
