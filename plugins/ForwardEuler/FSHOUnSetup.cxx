// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/FSHOUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FSHOUnSetup,
                      FwdEulerData,
                      ForwardEulerLib>
fshoUnSetupProvider("FSHOUnSetup");

//////////////////////////////////////////////////////////////////////////////

FSHOUnSetup::FSHOUnSetup(std::string name) : FwdEulerCom(name),
  socket_pastStates("pastStates"),
  socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

FSHOUnSetup::~FSHOUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FSHOUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastStates);
  result.push_back(&socket_interStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FSHOUnSetup::execute()
{
  CFAUTOTRACE;

  if (SubSystemStatusStack::getActive()->getDT() > 0.){
  CF_DEBUG_POINT;
    deleteAllPtr(socket_pastStates);
  CF_DEBUG_POINT;
    deleteAllPtr(socket_interStates);
  CF_DEBUG_POINT;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD
