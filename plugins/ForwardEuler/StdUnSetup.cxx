// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "ForwardEuler/ForwardEuler.hh"
#include "ForwardEuler/StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  

    namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup,
                      FwdEulerData,
                      ForwardEulerLib>
stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(std::string name) : FwdEulerCom(name),
  socket_pastStates("pastStates")
{
}

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::~StdUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdUnSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;

  if (SubSystemStatusStack::getActive()->getDT() > 0.)
  {
    deleteAllPtr(socket_pastStates);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ForwardEuler



} // namespace COOLFluiD
