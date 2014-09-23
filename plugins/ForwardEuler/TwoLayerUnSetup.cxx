// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ForwardEuler/ForwardEuler.hh"
#include "TwoLayerUnSetup.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerUnSetup, FwdEulerData, ForwardEulerLib> TwoLayerUnSetupProvider("TwoLayerUnSetup");

//////////////////////////////////////////////////////////////////////////////

TwoLayerUnSetup::TwoLayerUnSetup(std::string name) : FwdEulerCom(name),
socket_pastStates("pastStates"),
socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerUnSetup::~TwoLayerUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerUnSetup::execute()
{
  CFAUTOTRACE;

  DataHandle<State*> interStates = socket_interStates.getDataHandle();

  // deallocate the states in the storage
  for(CFuint i = 0; i < interStates.size(); ++i) {
    delete interStates[i];
  }

  if (SubSystemStatusStack::getActive()->getDT() > 0.){
    DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

    // deallocate the states in the storage
    for(CFuint i = 0; i < pastStates.size(); ++i) {
      delete pastStates[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD
